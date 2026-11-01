"""Resolve a public GSE/GSM record and cache its submitted processed H5MU.

NCBI's documented accession display endpoint supplies supplementary-file URLs.
Only NCBI HTTPS downloads are followed. No portal password or API key is needed
for released GEO data. Missing/private records and ambiguous files fail clearly.
"""
from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import re
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlparse, urlunparse, unquote
from urllib.request import HTTPRedirectHandler, Request, build_opener

from .igvf_download import digest_file

_RECORD = re.compile(r"^(GSE|GSM)[1-9][0-9]*$")
_HOSTS = {"www.ncbi.nlm.nih.gov", "ftp.ncbi.nlm.nih.gov"}
_EXTENSIONS = (".h5mu", ".h5mu.gz", ".h5ad", ".h5ad.gz")
_DATASET_TOKENS = {"hippo_timecourse": ("hippo",), "ipsc_k562_mix": ("k562",)}


def normalize_accession(accession):
    value = str(accession or "").strip().upper()
    if not _RECORD.fullmatch(value):
        raise ValueError("Enter a public GEO Series (GSE...) or Sample (GSM...) accession, not a URL or an SRA accession.")
    return value


def ncbi_url(url):
    parsed = urlparse(str(url).strip())
    if parsed.scheme not in {"https", "ftp"} or parsed.hostname not in _HOSTS:
        raise ValueError("GEO supplementary downloads must be official NCBI HTTPS/FTP URLs.")
    if parsed.username or parsed.password or parsed.port not in (None, 443) or not parsed.path.startswith("/geo/"):
        raise ValueError("Unexpected GEO supplementary URL.")
    return urlunparse(parsed._replace(scheme="https", fragment=""))


class NcbiRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return super().redirect_request(req, fp, code, msg, headers, ncbi_url(newurl))


class GeoClient:
    def __init__(self, *, opener=None, timeout=120):
        self.opener = opener or build_opener(NcbiRedirect())
        self.timeout = timeout

    def _open(self, url):
        url = ncbi_url(url)
        for attempt in range(3):
            try:
                return self.opener.open(Request(url, headers={"User-Agent": "cc-perturbseq-figure-reproduction"}), timeout=self.timeout)
            except HTTPError as exc:
                status = exc.code
                exc.close()
                if status in {401, 403, 404}:
                    raise FileNotFoundError("The GEO record/file is unavailable anonymously. Check the accession and release status, or use the local_file mode with your submitted H5MU. This workflow does not request your NCBI login password.") from None
                if status not in {429, 500, 502, 503, 504} or attempt == 2:
                    raise RuntimeError(f"GEO request failed with HTTP {status}.") from None
            except (URLError, TimeoutError):
                if attempt == 2:
                    raise RuntimeError("GEO could not be reached. Check network connectivity and retry; no partial download is accepted.") from None
            time.sleep(2 ** attempt)

    def _metadata(self, accession, target="self"):
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + urlencode({"acc": accession, "targ": target, "view": "brief", "form": "text"})
        with self._open(url) as response:
            raw = response.read(32 * 1024 * 1024 + 1)
        if len(raw) > 32 * 1024 * 1024:
            raise ValueError("GEO metadata exceeds the safety limit. Supply a specific GSM accession.")
        text = raw.decode("utf-8")
        if target == "self" and not re.search(r"^\^(?:SERIES|SAMPLE)\s*=\s*" + re.escape(accession) + r"\s*$", text, re.M):
            raise FileNotFoundError(f"GEO did not return a public SOFT record for {accession}. Check release status; use local_file mode before release.")
        result = []
        for line in text.splitlines():
            if re.match(r"^!(Series|Sample)_supplementary_file(?:_\d+)?\s*=", line):
                value = line.split("=", 1)[1].strip()
                if value.lower() not in {"none", "", "na"}:
                    result.append(ncbi_url(value))
        return sorted(set(result)), url

    def discover(self, accession, *, dataset, filename=None):
        accession = normalize_accession(accession)
        if dataset not in _DATASET_TOKENS:
            raise ValueError(f"Unknown figure dataset: {dataset}")
        urls, metadata_url = self._metadata(accession)
        candidates = [u for u in urls if urlparse(u).path.lower().endswith(_EXTENSIONS)]
        if not candidates and accession.startswith("GSE"):
            # Some submitters attach processed matrices to samples, not the series.
            sample_urls, _ = self._metadata(accession, target="gsm")
            candidates = [u for u in sample_urls if urlparse(u).path.lower().endswith(_EXTENSIONS)]
        if filename:
            if Path(str(filename)).name != str(filename):
                raise ValueError("GEO_FILENAME must be an exact supplementary filename, not a path.")
            candidates = [u for u in candidates if unquote(Path(urlparse(u).path).name) == filename]
        elif len(candidates) > 1:
            matched = [u for u in candidates if all(t in Path(urlparse(u).path).name.lower() for t in _DATASET_TOKENS[dataset])]
            if matched:
                candidates = matched
        if len(candidates) != 1:
            names = [unquote(Path(urlparse(u).path).name) for u in candidates]
            raise ValueError(f"Expected one processed matrix for {dataset} in {accession}; found {len(candidates)}: {names}. Set GEO_FILENAME to the exact H5MU filename or supply a specific GSM accession. Raw TAR/FASTQ archives are not downloaded as substitutes.")
        return {"accession": accession, "dataset": dataset, "metadata_url": metadata_url,
                "url": candidates[0], "filename": unquote(Path(urlparse(candidates[0]).path).name)}

    def download(self, accession, cache, *, dataset, filename=None, expected_sha256=None):
        import h5py
        record = self.discover(accession, dataset=dataset, filename=filename)
        name = record["filename"]
        if name != Path(name).name or name in {"", ".", ".."}:
            raise ValueError("Unsafe GEO filename.")
        if expected_sha256 is not None and not re.fullmatch(r"[a-fA-F0-9]{64}", str(expected_sha256)):
            raise ValueError("expected_sha256 must be the 64-character SHA-256 of the decompressed H5MU.")
        expected_sha256 = str(expected_sha256).lower() if expected_sha256 else None
        directory = Path(cache).expanduser().resolve() / record["accession"]
        directory.mkdir(parents=True, exist_ok=True)
        output = directory / (name[:-3] if name.lower().endswith(".gz") else name)
        receipt_path = output.with_name(output.name + ".download.json")
        if output.is_file() and receipt_path.is_file():
            receipt = json.loads(receipt_path.read_text())
            observed = digest_file(output, "sha256")
            if receipt.get("url") == record["url"] and observed == receipt.get("sha256") and (not expected_sha256 or observed == expected_sha256) and h5py.is_hdf5(output):
                return output
            # A corrupted or changed cache is never consumed.
        fd, temporary_name = tempfile.mkstemp(prefix=".geo-", dir=directory)
        os.close(fd)
        temporary = Path(temporary_name)
        plain = temporary.with_suffix(".hdf5")
        try:
            with self._open(record["url"]) as response, temporary.open("wb") as stream:
                expected_length = response.headers.get("Content-Length")
                copied = 0
                while block := response.read(4 * 1024 * 1024):
                    copied += len(block)
                    if copied > 64 * 1024 ** 3:
                        raise ValueError("GEO download exceeds 64 GiB safety limit.")
                    stream.write(block)
            if not copied or (expected_length is not None and copied != int(expected_length)):
                raise ValueError("Incomplete GEO download; file length does not match the response.")
            wire_sha256 = digest_file(temporary, "sha256")
            if name.lower().endswith(".gz"):
                with gzip.open(temporary, "rb") as source, plain.open("wb") as target:
                    size = 0
                    while block := source.read(4 * 1024 * 1024):
                        size += len(block)
                        if size > 64 * 1024 ** 3:
                            raise ValueError("Decompressed GEO file exceeds safety limit.")
                        target.write(block)
            else:
                temporary.rename(plain)
            if not h5py.is_hdf5(plain):
                raise ValueError("GEO returned neither an H5MU nor an H5AD container; no HTML or partial file is accepted.")
            observed = digest_file(plain, "sha256")
            if expected_sha256 and observed != expected_sha256:
                raise ValueError("GEO source SHA-256 mismatch; the downloaded object will not be used.")
            os.replace(plain, output)
            receipt = {**record, "sha256": observed, "download_sha256": wire_sha256,
                "size_bytes": output.stat().st_size, "expected_sha256_verified": bool(expected_sha256),
                "checksum_scope": "Caller-supplied source pin" if expected_sha256 else "Local download/cache integrity; no publisher checksum supplied"}
            receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
            return output
        finally:
            temporary.unlink(missing_ok=True)
            plain.unlink(missing_ok=True)

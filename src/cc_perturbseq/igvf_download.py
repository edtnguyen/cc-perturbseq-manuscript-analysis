"""Checksum-pinned IGVF downloads with optional, in-memory authentication."""
from __future__ import annotations

import base64
import getpass
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener

API = "https://api.data.igvf.org"


def digest_file(path, algorithm="md5"):
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class SafeRedirect(HTTPRedirectHandler):
    """Never forward IGVF credentials to the storage host in a redirect."""
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urlparse(newurl).scheme != "https":
            raise RuntimeError("Refusing a non-HTTPS download redirect.")
        redirected = super().redirect_request(req, fp, code, msg, headers, newurl)
        if redirected is not None and urlparse(newurl).netloc != urlparse(API).netloc:
            redirected.remove_header("Authorization")
        return redirected


class IgvfClient:
    def __init__(self, *, allow_prompt=True, timeout=120, opener=None):
        self.allow_prompt = allow_prompt
        self.timeout = timeout
        self.opener = opener or build_opener(SafeRedirect())
        self._key = os.environ.get("IGVF_ACCESS_KEY", "")
        self._secret = os.environ.get("IGVF_SECRET_ACCESS_KEY", "")
        if bool(self._key) != bool(self._secret):
            raise ValueError("Set both IGVF_ACCESS_KEY and IGVF_SECRET_ACCESS_KEY, or neither.")
        self._prompted = False

    def __repr__(self):
        return f"IgvfClient(auth_configured={bool(self._key)}, allow_prompt={self.allow_prompt})"

    def _open(self, url):
        parsed = urlparse(url)
        if parsed.scheme != "https" or parsed.netloc != urlparse(API).netloc:
            raise ValueError("Initial IGVF requests must use the official HTTPS API origin.")
        for attempt in range(4):
            request = Request(url, headers={"Accept": "application/json", "User-Agent": "cc-perturbseq-reproduction"})
            if self._key:
                encoded = base64.b64encode(f"{self._key}:{self._secret}".encode()).decode()
                request.add_header("Authorization", "Basic " + encoded)
            try:
                return self.opener.open(request, timeout=self.timeout)
            except HTTPError as exc:
                status = exc.code
                exc.close()
                if status in (401, 403):
                    if self.allow_prompt and not self._prompted:
                        self._prompted = True
                        print("IGVF denied anonymous/current access. Enter a portal access-key pair (not your login password).")
                        self._key = getpass.getpass("IGVF access key: ")
                        self._secret = getpass.getpass("IGVF secret access key: ")
                        if not self._key or not self._secret:
                            raise PermissionError("IGVF authentication was not provided.") from None
                        continue
                    raise PermissionError("IGVF access denied. Configure an authorized access-key pair; credentials are never saved.") from None
                if status in (429, 500, 502, 503, 504) and attempt < 3:
                    time.sleep(2 ** attempt)
                    continue
                raise RuntimeError(f"IGVF request failed (HTTP {status}); no input was substituted.") from None
            except URLError:
                if attempt < 3:
                    time.sleep(2 ** attempt)
                    continue
                raise RuntimeError("Could not reach the IGVF API. Check network access; no input was substituted.") from None
        raise RuntimeError("IGVF request exhausted retry limit.")

    def metadata(self, accession):
        with self._open(f"{API}/{accession}/?format=json") as response:
            return json.load(response)

    def download(self, record, cache_dir):
        accession = record["accession"]
        directory = Path(cache_dir).expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / record["filename"]
        receipt = directory / (record["filename"] + ".verified.json")
        expected = {"accession": accession, "md5sum": record["md5sum"], "file_size": int(record["file_size"])}
        if path.is_file() and path.stat().st_size == expected["file_size"]:
            stat = path.stat()
            current = {**expected, "mtime_ns": stat.st_mtime_ns}
            try:
                if json.loads(receipt.read_text()) == current:
                    return path
            except (OSError, ValueError):
                pass
            if digest_file(path) == expected["md5sum"]:
                receipt.write_text(json.dumps(current, indent=2) + "\n")
                return path
        metadata = self.metadata(accession)
        for key, value in expected.items():
            if metadata.get(key) != value:
                raise ValueError(f"IGVF metadata differs from the pinned {key} for {accession}; review the file version before use.")
        if metadata.get("upload_status") != "validated":
            raise ValueError(f"IGVF file {accession} is not upload-validated.")
        href = metadata.get("href")
        if not href:
            raise ValueError(f"IGVF file {accession} has no download endpoint.")
        url = urljoin(API, href)
        print(f"Downloading {accession} ({expected['file_size'] / 1024**2:.1f} MiB)")
        fd, temporary = tempfile.mkstemp(prefix=record["filename"] + ".", suffix=".part", dir=directory)
        md5 = hashlib.md5()
        size = 0
        try:
            with os.fdopen(fd, "wb") as target, self._open(url) as response:
                while True:
                    block = response.read(4 * 1024 * 1024)
                    if not block:
                        break
                    size += len(block)
                    if size > expected["file_size"]:
                        raise ValueError(f"Download exceeds pinned size: {accession}")
                    target.write(block)
                    md5.update(block)
            if size != expected["file_size"] or md5.hexdigest() != expected["md5sum"]:
                raise ValueError(f"Download failed size/MD5 verification: {accession}; cached input was not replaced.")
            os.replace(temporary, path)
            receipt.write_text(json.dumps({**expected, "mtime_ns": path.stat().st_mtime_ns}, indent=2) + "\n")
            return path
        finally:
            Path(temporary).unlink(missing_ok=True)

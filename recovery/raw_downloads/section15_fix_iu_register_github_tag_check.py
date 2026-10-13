IU_REGISTER_GITHUB_WRAPPER = (
    WORK_DIR / "iu_register_github_tags_wrapper.py"
)


def detect_igvf_utils_local_version() -> str:
    """Read the version from the external igvf_utils Python environment."""
    probe_code = r"""
from importlib.metadata import PackageNotFoundError, version

detected = ""
for distribution_name in ("igvf-utils", "igvf_utils"):
    try:
        detected = version(distribution_name)
        break
    except PackageNotFoundError:
        pass

if not detected:
    try:
        import igvf_utils
        detected = str(getattr(igvf_utils, "__version__", ""))
    except Exception:
        detected = ""

print(detected)
"""
    completed = subprocess.run(
        [str(IGVF_PYTHON), "-c", probe_code],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=os.environ.copy(),
        check=False,
    )

    version_candidates = [
        line.strip().lstrip("v")
        for line in completed.stdout.splitlines()
        if re.fullmatch(r"v?\d+(?:\.\d+)+", line.strip())
    ]

    if version_candidates:
        return version_candidates[-1]

    return "3.1.1"


def ensure_iu_register_github_wrapper() -> Path:
    """Create and syntax-check a narrow wrapper for the tags API call."""
    wrapper_source = '#!/usr/bin/env python\nfrom __future__ import annotations\n\nimport json\nimport os\nimport runpy\nimport sys\nfrom pathlib import Path\nfrom urllib.parse import urlparse\n\nimport requests\n\n\ndef _is_github_tags_url(url: object) -> bool:\n    try:\n        parsed = urlparse(str(url))\n    except Exception:\n        return False\n    hostname = (parsed.hostname or "").lower()\n    path = parsed.path.rstrip("/").lower()\n    return "github" in hostname and path.endswith("/tags")\n\n\nclass _FallbackTagResponse:\n    def __init__(self, local_version: str):\n        self._payload = [{"name": local_version}]\n        self.status_code = 200\n        self.ok = True\n        self.reason = "Local version fallback"\n        self.text = json.dumps(self._payload)\n        self.url = "local://igvf-utils-version-fallback"\n        self.headers = {}\n\n    def json(self):\n        return self._payload\n\n    def raise_for_status(self):\n        return None\n\n\ndef main() -> None:\n    if len(sys.argv) < 2:\n        raise SystemExit(\n            "Usage: iu_register_github_tags_wrapper.py "\n            "/path/to/iu_register.py [arguments...]"\n        )\n\n    target = Path(sys.argv[1]).resolve()\n    if not target.is_file():\n        raise FileNotFoundError(target)\n\n    local_version = os.environ.get(\n        "IGVF_UTILS_LOCAL_VERSION",\n        "3.1.1",\n    ).strip().lstrip("v")\n    if not local_version:\n        local_version = "3.1.1"\n\n    real_get = requests.get\n\n    def guarded_get(url, *args, **kwargs):\n        if not _is_github_tags_url(url):\n            return real_get(url, *args, **kwargs)\n\n        try:\n            response = real_get(url, *args, **kwargs)\n            payload = response.json()\n            valid = (\n                isinstance(payload, list)\n                and len(payload) > 0\n                and isinstance(payload[0], dict)\n                and bool(payload[0].get("name"))\n            )\n            if valid:\n                return response\n\n            preview = (\n                payload\n                if isinstance(payload, dict)\n                else str(payload)[:300]\n            )\n            print(\n                "[iu_register wrapper] GitHub tags response was not "\n                f"a usable tag list ({preview!r}); using installed "\n                f"local version {local_version} for the advisory "\n                "version check only.",\n                file=sys.stderr,\n            )\n        except Exception as error:\n            print(\n                "[iu_register wrapper] GitHub tags lookup failed "\n                f"({type(error).__name__}: {error}); using installed "\n                f"local version {local_version} for the advisory "\n                "version check only.",\n                file=sys.stderr,\n            )\n\n        return _FallbackTagResponse(local_version)\n\n    requests.get = guarded_get\n\n    sys.path.insert(0, str(target.parent))\n    sys.argv = [str(target), *sys.argv[2:]]\n    runpy.run_path(str(target), run_name="__main__")\n\n\nif __name__ == "__main__":\n    main()\n'

    if (
        not IU_REGISTER_GITHUB_WRAPPER.is_file()
        or IU_REGISTER_GITHUB_WRAPPER.read_text(encoding="utf-8")
        != wrapper_source
    ):
        IU_REGISTER_GITHUB_WRAPPER.write_text(
            wrapper_source,
            encoding="utf-8",
        )

    compile_result = subprocess.run(
        [
            str(IGVF_PYTHON),
            "-m",
            "py_compile",
            str(IU_REGISTER_GITHUB_WRAPPER),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=os.environ.copy(),
        check=False,
    )
    if compile_result.returncode != 0:
        raise RuntimeError(
            "The iu_register GitHub-tags wrapper failed syntax "
            "validation:\n"
            + compile_result.stdout
        )

    return IU_REGISTER_GITHUB_WRAPPER


def run_iu_register(
    *,
    profile: str,
    input_path: Path,
    dry_run: bool,
    patch: bool,
    label: str,
) -> dict:
    if not IGVF_PYTHON.is_file():
        raise FileNotFoundError(IGVF_PYTHON)
    if not IU_REGISTER.is_file():
        raise FileNotFoundError(IU_REGISTER)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if portal_auth() is None:
        raise RuntimeError("IGVF credentials are not set")

    wrapper_path = ensure_iu_register_github_wrapper()
    local_version = detect_igvf_utils_local_version()

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    mode = "dryrun" if dry_run else "submit"
    log_path = (
        LOG_DIR
        / f"{safe_stem(label)}_{mode}_{timestamp}.log"
    )

    command = [
        str(IGVF_PYTHON),
        str(wrapper_path),
        str(IU_REGISTER),
        "-m",
        "prod",
        "-p",
        profile,
        "--tries",
        "3",
        "--delay",
        "10",
        "--backoff",
        "2",
    ]

    if dry_run:
        command.append("-d")
    if patch:
        command.append("--patch")

    command.extend(["-i", str(input_path)])

    environment = os.environ.copy()
    environment["IGVF_UTILS_LOCAL_VERSION"] = local_version

    print("\nCommand:", " ".join(command))
    print("igvf_utils local version:", local_version)
    print("Log:", log_path)

    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
        check=False,
    )

    log_path.write_text(completed.stdout, encoding="utf-8")
    print(completed.stdout)

    fatal = re.search(
        (
            r"Traceback|ValidationFailure|"
            r"Failed to POST|Failed to PATCH|"
            r"\b(?:401|403|409|422) Client Error\b"
        ),
        completed.stdout,
        flags=re.IGNORECASE,
    )

    success = completed.returncode == 0 and fatal is None

    if not success:
        raise RuntimeError(
            "igvf_utils failed; review " f"{log_path}"
        )

    return {
        "success": True,
        "log_path": str(log_path),
        "igvf_utils_local_version": local_version,
        "used_github_tags_wrapper": True,
    }

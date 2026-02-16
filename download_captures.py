import argparse
import os
import re
import time
import urllib.error
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple


SERIES_MATRIX_DEFAULT = "GSE248728_series_matrix.txt"
DOWNLOAD_DIR_DEFAULT = "GSE248728_downloads"


@dataclass(frozen=True)
class RemoteFile:
    capture: str
    filetype: str
    url: str
    filename: str


def _ftp_to_https(url: str) -> str:
    # NCBI FTP is also reachable over HTTPS; HTTPS supports Range requests.
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return "https://ftp.ncbi.nlm.nih.gov/" + url[len("ftp://ftp.ncbi.nlm.nih.gov/") :]
    return url


def parse_series_matrix(series_matrix_path: str) -> Dict[str, Dict[str, RemoteFile]]:
    urls = []
    with open(series_matrix_path, "r") as f:
        for line in f:
            if line.startswith("!Series_supplementary_file"):
                match = re.search(r'"(ftp://[^"]+)"', line)
                if match:
                    urls.append(match.group(1))

    captures: Dict[str, Dict[str, RemoteFile]] = defaultdict(dict)
    for url in urls:
        url = _ftp_to_https(url)
        filename = os.path.basename(url)
        m = re.match(r"GSE248728_(capture\d+)_(barcodes|features|matrix)\.", filename)
        if not m:
            continue
        capture, filetype = m.group(1), m.group(2)
        captures[capture][filetype] = RemoteFile(
            capture=capture,
            filetype=filetype,
            url=url,
            filename=filename,
        )
    return captures


def _sorted_captures(captures: Iterable[str]) -> list[str]:
    return sorted(captures, key=lambda x: int(re.search(r"\d+", x).group()))


def _http_head(url: str, timeout: int, max_retries: int) -> Tuple[Optional[int], Optional[str]]:
    last_exc: Optional[BaseException] = None
    for attempt in range(max_retries + 1):
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                size = resp.headers.get("Content-Length")
                accept_ranges = resp.headers.get("Accept-Ranges")
                return (int(size) if size is not None else None, accept_ranges)
        except Exception as e:
            last_exc = e
            if attempt < max_retries:
                time.sleep(min(2**attempt, 30))
                continue
            break
    # HEAD failures are not fatal; we can still attempt a GET.
    return (None, None)


def _download_with_resume(
    url: str,
    dest_path: str,
    *,
    timeout: int,
    max_retries: int,
    expected_size: Optional[int],
) -> None:
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    part_path = dest_path + ".part"
    if os.path.exists(dest_path) and expected_size is not None and os.path.getsize(dest_path) == expected_size:
        return

    if os.path.exists(dest_path) and not os.path.exists(part_path):
        # If we don't know whether it's complete, treat it as a partial and resume safely.
        os.replace(dest_path, part_path)

    existing = os.path.getsize(part_path) if os.path.exists(part_path) else 0
    if expected_size is not None and existing > expected_size:
        # Local partial is larger than remote; start over.
        os.remove(part_path)
        existing = 0

    last_exc: Optional[BaseException] = None
    for attempt in range(max_retries + 1):
        try:
            headers = {}
            if existing > 0:
                headers["Range"] = f"bytes={existing}-"

            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                status = getattr(resp, "status", resp.getcode())
                if existing > 0 and status != 206:
                    # Server ignored Range; restart clean to avoid duplicating bytes.
                    existing = 0
                    if os.path.exists(part_path):
                        os.remove(part_path)

                mode = "ab" if existing > 0 else "wb"
                bytes_written = existing
                last_report = time.time()
                with open(part_path, mode) as out:
                    while True:
                        chunk = resp.read(1024 * 1024)
                        if not chunk:
                            break
                        out.write(chunk)
                        bytes_written += len(chunk)
                        now = time.time()
                        if now - last_report >= 5:
                            if expected_size:
                                pct = 100.0 * bytes_written / expected_size
                                print(f"    {os.path.basename(dest_path)}: {bytes_written}/{expected_size} bytes ({pct:.1f}%)")
                            else:
                                print(f"    {os.path.basename(dest_path)}: {bytes_written} bytes")
                            last_report = now

            final_size = os.path.getsize(part_path)
            if expected_size is not None and final_size != expected_size:
                raise IOError(f"Size mismatch after download: got {final_size}, expected {expected_size}")

            os.replace(part_path, dest_path)
            return
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ConnectionError, OSError) as e:
            last_exc = e
            if attempt < max_retries:
                wait = min(2**attempt, 30)
                print(f"    [retry {attempt + 1}/{max_retries}] {e} (waiting {wait}s)")
                time.sleep(wait)
                existing = os.path.getsize(part_path) if os.path.exists(part_path) else 0
                continue
            break
    raise RuntimeError(f"Failed to download after {max_retries} retries: {url}\nLast error: {last_exc}")


def _parse_csv_set(value: Optional[str]) -> Optional[set[str]]:
    if value is None:
        return None
    parts = [p.strip() for p in value.split(",") if p.strip()]
    return set(parts) if parts else None


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Download GSE248728 supplementary 10X matrices (barcodes/features/matrix) with resume + retries. "
            "Uses HTTPS mirror of NCBI GEO FTP to support Range requests."
        )
    )
    parser.add_argument("--series-matrix", default=SERIES_MATRIX_DEFAULT)
    parser.add_argument("--out-dir", default=DOWNLOAD_DIR_DEFAULT)
    parser.add_argument(
        "--captures",
        default=None,
        help='Comma-separated captures to download, e.g. "capture1,capture2". Default: all.',
    )
    parser.add_argument(
        "--types",
        default="barcodes,features,matrix",
        help='Comma-separated types: "barcodes,features,matrix" (default: all three).',
    )
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--retries", type=int, default=8)
    args = parser.parse_args()

    captures = parse_series_matrix(args.series_matrix)

    capture_filter = _parse_csv_set(args.captures)
    type_filter = _parse_csv_set(args.types) or {"barcodes", "features", "matrix"}

    selected_caps = _sorted_captures(
        cap for cap in captures.keys() if capture_filter is None or cap in capture_filter
    )
    if not selected_caps:
        raise SystemExit("No captures selected (check --captures)")

    print(f"Found {len(captures)} captures in series matrix; downloading {len(selected_caps)} capture(s).")
    os.makedirs(args.out_dir, exist_ok=True)
    manifest_path = os.path.join(args.out_dir, "manifest.tsv")

    with open(manifest_path, "w") as mf:
        mf.write("capture\tfiletype\turl\texpected_size\tlocal_size\tstatus\n")
        for cap in selected_caps:
            cap_dir = os.path.join(args.out_dir, cap)
            os.makedirs(cap_dir, exist_ok=True)
            files = captures[cap]
            for filetype in ["features", "barcodes", "matrix"]:
                if filetype not in type_filter:
                    continue
                remote = files.get(filetype)
                if remote is None:
                    print(f"  WARNING: {cap} missing {filetype}")
                    mf.write(f"{cap}\t{filetype}\t\t\t\tmissing\n")
                    continue

                dest = os.path.join(cap_dir, remote.filename)
                expected_size, _ = _http_head(remote.url, timeout=args.timeout, max_retries=max(1, args.retries // 2))

                try:
                    print(f"  {cap}: {remote.filename}")
                    _download_with_resume(
                        remote.url,
                        dest,
                        timeout=args.timeout,
                        max_retries=args.retries,
                        expected_size=expected_size,
                    )
                    local_size = os.path.getsize(dest)
                    mf.write(
                        f"{cap}\t{filetype}\t{remote.url}\t{expected_size if expected_size is not None else ''}\t{local_size}\tok\n"
                    )
                except Exception as e:
                    local_size = os.path.getsize(dest + ".part") if os.path.exists(dest + ".part") else ""
                    error_msg = str(e).replace("\t", " ")
                    mf.write(
                        f"{cap}\t{filetype}\t{remote.url}\t{expected_size if expected_size is not None else ''}\t{local_size}\terror: {error_msg}\n"
                    )
                    print(f"    [ERROR] {e}")

    print(f"\nDone. Wrote manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

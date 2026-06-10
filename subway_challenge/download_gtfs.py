"""Download and extract the MTA NYC Subway static GTFS feed.

The canonical static feed lives at:
    https://rrgtfsfeeds.s3.amazonaws.com/gtfs_subway.zip
(the legacy ``web.mta.info/developers/data/nyct/subway/google_transit.zip``
URL 301-redirects to the same S3 object).

Usage:
    python -m subway_challenge.download_gtfs            # download + extract
    python -m subway_challenge.download_gtfs --force    # re-download
"""
from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path

import requests

GTFS_SUBWAY_URL = "https://rrgtfsfeeds.s3.amazonaws.com/gtfs_subway.zip"

# Files we expect in a valid subway GTFS feed.
REQUIRED_FILES = {"stops.txt", "stop_times.txt", "trips.txt", "routes.txt", "calendar.txt"}

DEFAULT_DEST = Path("data/gtfs")


def download_gtfs(
    dest_dir: Path | str = DEFAULT_DEST,
    url: str = GTFS_SUBWAY_URL,
    force: bool = False,
    chunk_size: int = 1 << 16,
) -> Path:
    """Download the GTFS zip and extract it into ``dest_dir``.

    Returns the directory the feed was extracted into. Skips the download if the
    feed already appears present unless ``force`` is set.
    """
    dest_dir = Path(dest_dir)

    if not force and _feed_present(dest_dir):
        print(f"GTFS already present at {dest_dir} (use force=True to re-download).")
        return dest_dir

    dest_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading GTFS from {url} ...")

    resp = requests.get(url, stream=True, timeout=60)
    resp.raise_for_status()

    total = int(resp.headers.get("Content-Length", 0))
    buf = io.BytesIO()
    downloaded = 0
    for chunk in resp.iter_content(chunk_size=chunk_size):
        buf.write(chunk)
        downloaded += len(chunk)
        if total:
            pct = 100 * downloaded / total
            print(f"\r  {downloaded/1e6:5.1f} / {total/1e6:5.1f} MB ({pct:4.1f}%)", end="")
    print()

    buf.seek(0)
    with zipfile.ZipFile(buf) as zf:
        names = set(zf.namelist())
        missing = REQUIRED_FILES - names
        if missing:
            raise RuntimeError(
                f"Downloaded archive is missing expected GTFS files: {sorted(missing)}"
            )
        zf.extractall(dest_dir)

    last_modified = resp.headers.get("Last-Modified", "unknown")
    print(f"Extracted {len(names)} files to {dest_dir} (feed Last-Modified: {last_modified}).")
    return dest_dir


def _feed_present(dest_dir: Path) -> bool:
    return dest_dir.is_dir() and all((dest_dir / f).exists() for f in REQUIRED_FILES)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Download MTA NYC Subway GTFS feed.")
    parser.add_argument("--dest", default=str(DEFAULT_DEST), help="Destination directory.")
    parser.add_argument("--url", default=GTFS_SUBWAY_URL, help="GTFS zip URL.")
    parser.add_argument("--force", action="store_true", help="Re-download even if present.")
    args = parser.parse_args(argv)

    download_gtfs(dest_dir=args.dest, url=args.url, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())

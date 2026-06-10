"""Canonical station identity for the Subway Challenge.

GTFS splits a few physical stations into multiple parent stations (e.g. the
upper/lower levels of W 4 St). The MTA "Subway Stations" open dataset assigns
each GTFS parent an official ``Station ID``; collapsing on it yields the
**472** stations the Subway Challenge counts (475 GTFS parents - 3 bi-level
splits, after dropping the 21 Staten Island Railway stops).

This module maps any GTFS stop (platform ``A12N`` or parent ``A12``) to its
official Station ID, and exposes the canonical 472-station set so a solver can
check "visited every station" exactly.

The mapping uses three columns of the inventory:
* ``GTFS Stop ID``  -- the GTFS parent station id (our ``station_of`` output)
* ``Station ID``    -- one physical station (the 472 count)  <-- canonical here
* ``Complex ID``    -- a fare-controlled complex (the 424 count; NOT used as the
                       station identity -- it over-merges Times Sq etc.)

Download the inventory with :func:`download_official_stations` (or it is fetched
on demand by :meth:`StationIndex.load`).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from .build_graph import station_of

OFFICIAL_URL = "https://data.ny.gov/api/views/39hk-dx4f/rows.csv?accessType=DOWNLOAD"
DEFAULT_PATH = Path("data/official/mta_subway_stations.csv")
SIR_DIVISION = "SIR"


def download_official_stations(path: Path | str = DEFAULT_PATH) -> Path:
    """Download the MTA Subway Stations inventory CSV to ``path``."""
    import requests

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(OFFICIAL_URL, timeout=60)
    resp.raise_for_status()
    path.write_bytes(resp.content)
    print(f"Downloaded MTA station inventory -> {path} ({len(resp.content)} bytes)")
    return path


class StationIndex:
    """Maps GTFS stops to official MTA Station IDs (subway only, SIR excluded)."""

    def __init__(self, inventory: pd.DataFrame, include_sir: bool = False):
        df = inventory if include_sir else inventory[inventory["Division"] != SIR_DIVISION]

        # GTFS parent station id -> official Station ID (and friendly metadata)
        self.parent_to_station: dict[str, str] = dict(
            zip(df["GTFS Stop ID"], df["Station ID"])
        )
        self.station_name: dict[str, str] = (
            df.groupby("Station ID")["Stop Name"].first().to_dict()
        )
        self.station_complex: dict[str, str] = (
            df.groupby("Station ID")["Complex ID"].first().to_dict()
        )
        self.station_to_parents: dict[str, list[str]] = (
            df.groupby("Station ID")["GTFS Stop ID"].apply(list).to_dict()
        )
        # The canonical universe: every distinct official Station ID.
        self.canonical_stations: frozenset[str] = frozenset(df["Station ID"])

    @classmethod
    def load(cls, path: Path | str = DEFAULT_PATH, include_sir: bool = False,
             download_if_missing: bool = True) -> "StationIndex":
        path = Path(path)
        if not path.exists():
            if not download_if_missing:
                raise FileNotFoundError(f"{path} not found; call download_official_stations().")
            download_official_stations(path)
        return cls(pd.read_csv(path, dtype=str), include_sir=include_sir)

    def resolve(self, stop_id: str) -> str:
        """Official Station ID for a platform/parent/GTFS stop id.

        Falls back to the GTFS parent id itself if the stop is absent from the
        inventory (keeps it a distinct station rather than dropping it).
        """
        parent = station_of(stop_id)
        return self.parent_to_station.get(parent, parent)

    def collapses(self) -> dict[str, list[str]]:
        """Official stations that map to >1 GTFS parent (the bi-level splits)."""
        return {sid: ps for sid, ps in self.station_to_parents.items() if len(ps) > 1}

    def __len__(self) -> int:
        return len(self.canonical_stations)


def main(argv: list[str] | None = None) -> int:
    idx = StationIndex.load()
    print(f"Canonical subway stations (official Station IDs, SIR excluded): {len(idx)}")
    print(f"Distinct GTFS parents mapped: {len(idx.parent_to_station)}")
    print(f"Distinct complexes: {len(set(idx.station_complex.values()))}")
    print("\nBi-level stations collapsed (one official station = two GTFS parents):")
    for sid, parents in idx.collapses().items():
        print(f"  Station {sid} '{idx.station_name[sid]}': {parents}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Out-of-system run distances between station complexes, via hosted OSRM.

Provides the run/walk data the on-demand :class:`~subway_challenge.run_layer.RunLayer`
is built from: complex-to-complex street distances (cached once from the public
OSRM demo server) converted to run times via a pace + per-end access penalty.

`load_complexes` reads the MTA station + entrance datasets; `osrm_full_matrix`
fetches/caches the all-pairs street-distance matrix; `complex_run_adjacency`
turns it into `{complex: [(other_complex, seconds, meters), ...]}`, excluding
same-complex and in-system-linked pairs.
"""
from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

DEFAULT_OSRM_HOST = "http://router.project-osrm.org"
DEFAULT_OSRM_PROFILE = "driving"          # public demo only serves driving; we use distance
# Run model. Lenient settings (committed): consistent with real record-holders,
# whose aggressive running puts the record (22:14) above this model's lower bound
# (22:01). A stricter model (2.5 m/s, 1.5 km cap) was rejected: its lower bound
# (22:22) exceeded the record -- i.e. provably too strict.
DEFAULT_PACE_MPS = 3.0                     # ~10.8 km/h
DEFAULT_ACCESS_PENALTY_S = 60             # per-end station ingress/egress
DEFAULT_RUN_MAX_METERS = None             # no single-run distance cap
DEFAULT_RUN_RADIUS_M = 5000.0

DEFAULT_GTFS = Path("data/gtfs")
STATIONS_CSV = Path("data/official/mta_subway_stations.csv")
ENTRANCES_CSV = Path("data/official/mta_subway_entrances.csv")
OSRM_CACHE = Path("data/walk/osrm_matrix.json")


@dataclass
class Complex:
    cid: str
    parents: list[str]
    lat: float                  # centroid
    lon: float
    exits: list[tuple]          # (lat, lon) entrances with Exit Allowed
    entries: list[tuple]        # (lat, lon) entrances with Entry Allowed


def _haversine_m(lat1, lon1, lat2, lon2) -> float:
    r = 6_371_000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def load_complexes(stations_csv=STATIONS_CSV, entrances_csv=ENTRANCES_CSV,
                   include_sir: bool = False) -> dict[str, Complex]:
    """Load station complexes (subway only) with their parents and street
    entrances (exit-allowed / entry-allowed), keyed by official Complex ID."""
    st = pd.read_csv(stations_csv, dtype=str)
    en = pd.read_csv(entrances_csv, dtype=str)
    if not include_sir:
        st = st[st["Division"] != "SIR"]
        en = en[en["Division"] != "SIR"]
    st = st.copy()
    st["lat"] = st["GTFS Latitude"].astype(float)
    st["lon"] = st["GTFS Longitude"].astype(float)
    en = en.copy()
    en["lat"] = en["Entrance Latitude"].astype(float)
    en["lon"] = en["Entrance Longitude"].astype(float)

    def _points_by_complex(df) -> dict[str, list]:
        out: dict[str, list] = defaultdict(list)
        for cid, lat, lon in zip(df["Complex ID"], df["lat"], df["lon"]):
            out[cid].append((lat, lon))
        return out

    ent_by_cx_exit = _points_by_complex(en[en["Exit Allowed"].str.upper() != "NO"])
    ent_by_cx_entry = _points_by_complex(en[en["Entry Allowed"].str.upper() != "NO"])

    complexes: dict[str, Complex] = {}
    for cid, g in st.groupby("Complex ID"):
        clat, clon = g["lat"].mean(), g["lon"].mean()
        exits = ent_by_cx_exit.get(cid) or [(clat, clon)]      # fall back to centroid
        entries = ent_by_cx_entry.get(cid) or [(clat, clon)]
        complexes[cid] = Complex(cid, list(g["GTFS Stop ID"]), clat, clon, exits, entries)
    return complexes


def _osrm_table(coords, src_idx, dst_idx, host, profile, session, timeout=60):
    """One OSRM /table call. Returns a len(src_idx) x len(dst_idx) distance grid
    (meters). ``coords`` are (lat, lon); src_idx/dst_idx are global indices."""
    uniq = list(dict.fromkeys(src_idx + dst_idx))      # dedupe, keep order
    pos = {g: p for p, g in enumerate(uniq)}
    coordstr = ";".join(f"{coords[g][1]},{coords[g][0]}" for g in uniq)  # lon,lat
    sources = ";".join(str(pos[g]) for g in src_idx)
    dests = ";".join(str(pos[g]) for g in dst_idx)
    url = (f"{host}/table/v1/{profile}/{coordstr}"
           f"?annotations=distance&sources={sources}&destinations={dests}")
    r = session.get(url, timeout=timeout)
    r.raise_for_status()
    data = r.json()
    if data.get("code") != "Ok":
        raise RuntimeError(f"OSRM error: {data.get('code')} {data.get('message','')}")
    return data["distances"]


def osrm_full_matrix(coords, host=DEFAULT_OSRM_HOST, profile=DEFAULT_OSRM_PROFILE,
                     block=50, cache_path=OSRM_CACHE, delay=0.3):
    """Full NxN street-distance matrix (meters) over ``coords`` via batched OSRM
    /table calls, cached to disk so re-runs are free. ``None`` = no route."""
    import time
    import requests

    n = len(coords)
    cache = json.loads(Path(cache_path).read_text()) if Path(cache_path).exists() else {}

    def key(i, j):
        return f"{coords[i][0]:.6f},{coords[i][1]:.6f}|{coords[j][0]:.6f},{coords[j][1]:.6f}"

    D = [[0.0 if i == j else None for j in range(n)] for i in range(n)]
    blocks = [list(range(s, min(s + block, n))) for s in range(0, n, block)]
    session = requests.Session()
    calls = 0
    for bi in blocks:
        for bj in blocks:
            need = [(i, j) for i in bi for j in bj if i != j and key(i, j) not in cache]
            if need:
                grid = _osrm_table(coords, bi, bj, host, profile, session)
                for a, i in enumerate(bi):
                    for b, j in enumerate(bj):
                        if i != j:
                            cache[key(i, j)] = grid[a][b]
                calls += 1
                if delay:
                    time.sleep(delay)
            for i in bi:
                for j in bj:
                    if i != j:
                        D[i][j] = cache.get(key(i, j))
    Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
    Path(cache_path).write_text(json.dumps(cache))
    print(f"OSRM matrix: {n} points, {calls} new table calls (cache now {len(cache)} cells).")
    return D


def complex_run_adjacency(radius_m=DEFAULT_RUN_RADIUS_M, pace_mps=DEFAULT_PACE_MPS,
                          access_penalty_s=DEFAULT_ACCESS_PENALTY_S, min_seconds=0,
                          max_seconds=None, max_meters=DEFAULT_RUN_MAX_METERS,
                          host=DEFAULT_OSRM_HOST, profile=DEFAULT_OSRM_PROFILE,
                          gtfs_dir=DEFAULT_GTFS):
    """Complex-to-complex run adjacency within ``radius_m``: returns
    ``(adjacency, complexes)`` where adjacency[cid] = [(other_cid, seconds, meters), ...].

    Street distance (symmetrized across one-way asymmetry, from the cached OSRM
    matrix) -> time via ``pace_mps`` + per-end ``access_penalty_s``. Excludes
    same-complex and pairs already linked in-system by transfers.txt.
    """
    complexes = load_complexes()
    cids = sorted(complexes)
    coords = [(complexes[c].lat, complexes[c].lon) for c in cids]

    tr = pd.read_csv(Path(gtfs_dir) / "transfers.txt", dtype=str)
    linked = {frozenset((a, b)) for a, b in zip(tr["from_stop_id"], tr["to_stop_id"]) if a != b}

    D = osrm_full_matrix(coords, host=host, profile=profile)

    adjacency: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    n = len(cids)
    for i in range(n):
        ca = complexes[cids[i]]
        for j in range(i + 1, n):
            cb = complexes[cids[j]]
            if _haversine_m(ca.lat, ca.lon, cb.lat, cb.lon) > radius_m:
                continue
            if any(frozenset((pa, pb)) in linked for pa in ca.parents for pb in cb.parents):
                continue
            cands = [d for d in (D[i][j], D[j][i]) if d is not None]
            if not cands:
                continue
            meters = int(round(min(cands)))
            if max_meters is not None and meters > max_meters:
                continue
            secs = max(int(round(min(cands) / pace_mps)) + access_penalty_s, min_seconds)
            if max_seconds is not None and secs > max_seconds:
                continue
            adjacency[cids[i]].append((cids[j], secs, meters))
            adjacency[cids[j]].append((cids[i], secs, meters))
    return adjacency, complexes

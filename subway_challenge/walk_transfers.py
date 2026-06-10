"""Precise street-level walking transfers between nearby, unlinked complexes.

GTFS ``transfers.txt`` covers in-system transfers (within a fare-controlled
complex). This module is the separate, out-of-system piece: walking links
between two *different* station complexes that sit a short walk apart but share
no free transfer (e.g. a block between complexes a challenger would walk).

Design (matches the chosen "Google walking, entrance-to-entrance" approach):

* Work at the **complex** level. All of a complex's street entrances serve all
  its lines, and an out-of-system walk is really complex-to-complex. This also
  fixes the 14 GTFS parents whose entrances are filed under a sibling parent.
* Endpoints are real **street entrances** (MTA "Subway Entrances and Exits"
  dataset), not platform centroids. For a walk CA -> CB we depart an
  *exit-allowed* entrance of CA and arrive an *entry-allowed* entrance of CB.
* For each directed complex pair we route the ``k`` geometrically-closest
  entrance pairs and keep the **minimum** walking duration.
* The engine is Google's **Routes API** (``computeRouteMatrix``, ``WALK``).
  Every element is cached to disk so each pair is billed at most once.
* Results are exported at **GTFS-parent granularity** (every parent of CA to
  every parent of CB) so ``build_graph`` can add them exactly like
  ``transfers.txt`` links, as ``mode="transfer", walk=True`` edges.

Cost: ~$5 / 1000 elements (Routes Essentials). At radius 400 m, k=3 that is
~2.1k elements (~$11) one-off, then cached. Use ``--estimate`` to price first.

ToS note: Google restricts long-term storage of results. The on-disk cache is a
local compute cache; treat it accordingly. (OSRM/Valhalla self-host avoid this.)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

WALK_SPEED_MPS = 1.35           # used only by the offline "dummy" backend
DEFAULT_RADIUS_M = 400.0
DEFAULT_K = 3                   # nearest entrance pairs routed per directed pair
GOOGLE_PRICE_PER_ELEMENT = 0.005

# OSRM (hosted) run-transfer defaults. The public demo only serves the driving
# profile, so we take its *distance* (street-network meters, a good proxy for a
# runner's path) and convert to time with our own pace -- never trusting the
# profile's speed.
DEFAULT_OSRM_HOST = "http://router.project-osrm.org"
DEFAULT_OSRM_PROFILE = "driving"
DEFAULT_PACE_MPS = 3.0          # ~10.8 km/h running pace
DEFAULT_ACCESS_PENALTY_S = 60   # per-end station ingress/egress not in street time
DEFAULT_RUN_RADIUS_M = 5000.0

DEFAULT_GTFS = Path("data/gtfs")
STATIONS_CSV = Path("data/official/mta_subway_stations.csv")
ENTRANCES_CSV = Path("data/official/mta_subway_entrances.csv")
OUT_CSV = Path("data/walk/walk_transfers.csv")
CACHE_JSON = Path("data/walk/google_cache.json")
OSRM_CACHE = Path("data/walk/osrm_matrix.json")

ROUTES_MATRIX_URL = "https://routes.googleapis.com/distanceMatrix/v2:computeRouteMatrix"


@dataclass(frozen=True)
class WalkLink:
    from_parent: str
    to_parent: str
    seconds: int
    meters: int
    from_complex: str
    to_complex: str


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


# -- load complexes ----------------------------------------------------------

def load_complexes(stations_csv=STATIONS_CSV, entrances_csv=ENTRANCES_CSV,
                   include_sir: bool = False) -> dict[str, Complex]:
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


def candidate_complex_pairs(complexes, radius_m=DEFAULT_RADIUS_M,
                            gtfs_dir=DEFAULT_GTFS):
    """Unordered complex pairs whose centroids are within ``radius_m`` and that
    are not already linked in-system by ``transfers.txt``."""
    tr = pd.read_csv(Path(gtfs_dir) / "transfers.txt", dtype=str)
    linked = {frozenset((a, b)) for a, b in zip(tr["from_stop_id"], tr["to_stop_id"]) if a != b}
    parent_cx = {p: c.cid for c in complexes.values() for p in c.parents}

    def is_linked(ca, cb):
        return any(frozenset((pa, pb)) in linked
                   for pa in ca.parents for pb in cb.parents)

    items = list(complexes.values())
    out = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            ca, cb = items[i], items[j]
            d = _haversine_m(ca.lat, ca.lon, cb.lat, cb.lon)
            if d <= radius_m and not is_linked(ca, cb):
                out.append((ca, cb, d))
    return out


def _nearest(points, toward, k):
    return sorted(points, key=lambda p: _haversine_m(p[0], p[1], toward[0], toward[1]))[:k]


# -- Google Routes backend (with disk cache) ---------------------------------

class GoogleWalk:
    def __init__(self, api_key: str, cache_path: Path = CACHE_JSON):
        import requests
        self._requests = requests
        self.api_key = api_key
        self.cache_path = Path(cache_path)
        self.cache: dict[str, dict] = {}
        if self.cache_path.exists():
            self.cache = json.loads(self.cache_path.read_text())
        self.billed_elements = 0

    @staticmethod
    def _key(o, d) -> str:
        return f"{o[0]:.6f},{o[1]:.6f}|{d[0]:.6f},{d[1]:.6f}"

    def matrix(self, origins, destinations):
        """Return {(oi, di): (seconds, meters)} for routable elements, using and
        updating the cache; only un-cached elements are sent to Google."""
        result, missing = {}, []
        for oi, o in enumerate(origins):
            for di, d in enumerate(destinations):
                c = self.cache.get(self._key(o, d))
                if c is not None:
                    if c.get("ok"):
                        result[(oi, di)] = (c["s"], c["m"])
                else:
                    missing.append((oi, di))
        if missing:
            self._query(origins, destinations, missing, result)
        return result

    def _query(self, origins, destinations, missing, result):
        # Routes matrix is rectangular; request the full grid covering missing
        # elements (simplest correct call), then cache every returned element.
        oset = sorted({oi for oi, _ in missing})
        dset = sorted({di for _, di in missing})
        body = {
            "origins": [{"waypoint": {"location": {"latLng":
                        {"latitude": origins[oi][0], "longitude": origins[oi][1]}}}}
                        for oi in oset],
            "destinations": [{"waypoint": {"location": {"latLng":
                            {"latitude": destinations[di][0], "longitude": destinations[di][1]}}}}
                            for di in dset],
            "travelMode": "WALK",
        }
        headers = {
            "Content-Type": "application/json",
            "X-Goog-Api-Key": self.api_key,
            "X-Goog-FieldMask": "originIndex,destinationIndex,duration,distanceMeters,condition",
        }
        resp = self._requests.post(ROUTES_MATRIX_URL, json=body, headers=headers, timeout=60)
        resp.raise_for_status()
        for el in resp.json():
            oi, di = oset[el.get("originIndex", 0)], dset[el.get("destinationIndex", 0)]
            o, d = origins[oi], destinations[di]
            self.billed_elements += 1
            if el.get("condition") == "ROUTE_EXISTS" and "duration" in el:
                s = int(str(el["duration"]).rstrip("s"))
                m = int(el.get("distanceMeters", 0))
                self.cache[self._key(o, d)] = {"ok": True, "s": s, "m": m}
                result[(oi, di)] = (s, m)
            else:
                self.cache[self._key(o, d)] = {"ok": False}

    def save(self):
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.cache_path.write_text(json.dumps(self.cache))


def _directed_walk(backend, src: Complex, dst: Complex, k: int):
    """Min walking (seconds, meters) from an exit of ``src`` to an entry of
    ``dst`` over the k nearest entrance pairs. ``None`` if no route exists."""
    origins = _nearest(src.exits, (dst.lat, dst.lon), k)
    dests = _nearest(dst.entries, (src.lat, src.lon), k)
    if backend == "dummy":
        best = min(((_haversine_m(*o, *d), o, d) for o in origins for d in dests),
                   key=lambda t: t[0])
        return int(round(best[0] / WALK_SPEED_MPS)), int(round(best[0]))
    cells = backend.matrix(origins, dests)        # GoogleWalk instance
    if not cells:
        return None
    s, m = min(cells.values(), key=lambda t: t[0])
    return s, m


def compute_walk_links(backend_name="google", radius_m=DEFAULT_RADIUS_M, k=DEFAULT_K,
                       max_seconds=None, gtfs_dir=DEFAULT_GTFS, out_csv=OUT_CSV,
                       api_key=None) -> list[WalkLink]:
    complexes = load_complexes()
    pairs = candidate_complex_pairs(complexes, radius_m, gtfs_dir)

    if backend_name == "google":
        key = api_key or os.environ.get("GOOGLE_MAPS_API_KEY")
        if not key:
            raise RuntimeError("Set GOOGLE_MAPS_API_KEY (or pass api_key=) for the google backend.")
        backend = GoogleWalk(key)
    elif backend_name == "dummy":
        backend = "dummy"
    else:
        raise ValueError(f"unknown backend {backend_name!r}")

    links: list[WalkLink] = []
    for ca, cb, _ in pairs:
        for src, dst in ((ca, cb), (cb, ca)):
            res = _directed_walk(backend, src, dst, k)
            if res is None:
                continue
            secs, meters = res
            if max_seconds is not None and secs > max_seconds:
                continue
            for pa in src.parents:
                for pb in dst.parents:
                    links.append(WalkLink(pa, pb, secs, meters, src.cid, dst.cid))

    if isinstance(backend, GoogleWalk):
        backend.save()
        print(f"Google billed {backend.billed_elements} new elements "
              f"(~${backend.billed_elements * GOOGLE_PRICE_PER_ELEMENT:.2f}); "
              f"cache now {len(backend.cache)} elements.")

    _write_csv(links, out_csv)
    return links


def _write_csv(links, out_csv):
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([l.__dict__ for l in links]).to_csv(out_csv, index=False)
    print(f"Wrote {len(links)} directed parent-level walk links -> {out_csv}")


# -- OSRM (hosted) full-matrix run transfers ---------------------------------

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
            if not need:
                pass
            else:
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
    print(f"OSRM matrix: {n} points, {calls} new table calls "
          f"(cache now {len(cache)} cells).")
    return D


def complex_run_adjacency(radius_m=DEFAULT_RUN_RADIUS_M, pace_mps=DEFAULT_PACE_MPS,
                          access_penalty_s=DEFAULT_ACCESS_PENALTY_S, min_seconds=0,
                          max_seconds=None, host=DEFAULT_OSRM_HOST,
                          profile=DEFAULT_OSRM_PROFILE, gtfs_dir=DEFAULT_GTFS):
    """Complex-to-complex run adjacency within ``radius_m``: returns
    ``(adjacency, complexes)`` where adjacency[cid] = [(other_cid, seconds), ...].

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

    adjacency: dict[str, list[tuple[str, int]]] = defaultdict(list)
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
            secs = max(int(round(min(cands) / pace_mps)) + access_penalty_s, min_seconds)
            if max_seconds is not None and secs > max_seconds:
                continue
            adjacency[cids[i]].append((cids[j], secs, meters))
            adjacency[cids[j]].append((cids[i], secs, meters))
    return adjacency, complexes


def compute_run_links(radius_m=DEFAULT_RUN_RADIUS_M, pace_mps=DEFAULT_PACE_MPS,
                      access_penalty_s=DEFAULT_ACCESS_PENALTY_S, min_seconds=0,
                      max_seconds=None, host=DEFAULT_OSRM_HOST,
                      profile=DEFAULT_OSRM_PROFILE, gtfs_dir=DEFAULT_GTFS,
                      out_csv=OUT_CSV) -> list[WalkLink]:
    """Export parent-level run links to CSV (for baking into a static graph).
    For on-demand use in the solver, prefer :func:`complex_run_adjacency`."""
    adjacency, complexes = complex_run_adjacency(
        radius_m, pace_mps, access_penalty_s, min_seconds, max_seconds,
        host, profile, gtfs_dir)
    seen, links = set(), []
    for cid, nbrs in adjacency.items():
        ca = complexes[cid]
        for ocid, secs, meters in nbrs:
            if (cid, ocid) in seen:
                continue
            seen.add((cid, ocid))
            cb = complexes[ocid]
            for pa in ca.parents:
                for pb in cb.parents:
                    links.append(WalkLink(pa, pb, secs, meters, cid, ocid))
    _write_csv(links, out_csv)
    print(f"radius={radius_m:.0f}m pace={pace_mps} m/s penalty={access_penalty_s}s "
          f"-> {len(links)} directed run links")
    return links


def estimate(radius_m=DEFAULT_RADIUS_M, k=DEFAULT_K):
    complexes = load_complexes()
    pairs = candidate_complex_pairs(complexes, radius_m)
    elems = 0
    for ca, cb, _ in pairs:
        for src, dst in ((ca, cb), (cb, ca)):
            elems += min(k, len(src.exits)) * min(k, len(dst.entries))
    print(f"radius={radius_m:.0f}m k={k}: {len(pairs)} complex pairs, "
          f"~{elems} Google elements -> ~${elems * GOOGLE_PRICE_PER_ELEMENT:.2f} "
          f"(worst case, before cache hits)")
    return elems


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Compute walk/run transfer links.")
    p.add_argument("--backend", default="osrm", choices=["osrm", "google", "dummy"],
                   help="osrm: hosted OSRM run links (free, default). "
                        "google: precise entrance-to-entrance walks (paid).")
    p.add_argument("--radius", type=float, default=None,
                   help="Max straight-line distance to link (default 5000 for osrm, 400 for google).")
    p.add_argument("--k", type=int, default=DEFAULT_K, help="(google) nearest entrance pairs.")
    p.add_argument("--pace", type=float, default=DEFAULT_PACE_MPS,
                   help="(osrm) running pace in m/s (default 3.0 ~ 10.8 km/h).")
    p.add_argument("--access-penalty", type=int, default=DEFAULT_ACCESS_PENALTY_S,
                   help="(osrm) per-end station ingress/egress seconds (default 60).")
    p.add_argument("--max-seconds", type=int, default=None,
                   help="Drop links whose transfer exceeds this many seconds.")
    p.add_argument("--osrm-host", default=DEFAULT_OSRM_HOST)
    p.add_argument("--profile", default=DEFAULT_OSRM_PROFILE)
    p.add_argument("--out", default=str(OUT_CSV))
    p.add_argument("--estimate", action="store_true",
                   help="(google) print cost estimate and exit (no API calls).")
    args = p.parse_args(argv)

    if args.estimate:
        estimate(args.radius or DEFAULT_RADIUS_M, args.k)
        return 0
    if args.backend == "osrm":
        compute_run_links(radius_m=args.radius or DEFAULT_RUN_RADIUS_M, pace_mps=args.pace,
                          access_penalty_s=args.access_penalty, max_seconds=args.max_seconds,
                          host=args.osrm_host, profile=args.profile, out_csv=args.out)
    else:
        compute_walk_links(args.backend, args.radius or DEFAULT_RADIUS_M, args.k,
                           args.max_seconds, out_csv=args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

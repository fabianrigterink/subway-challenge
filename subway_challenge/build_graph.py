"""Build a time-expanded graph of the NYC subway from GTFS.

Model
-----
A *node* is an event ``(stop_id, t)`` where ``stop_id`` is a platform-level GTFS
stop (e.g. ``127N``) and ``t`` is a time in seconds within a cyclic week,
``t in [0, 604800)`` with Monday 00:00:00 == 0. After-midnight GTFS times
(``>= 24:00:00``) and Sunday-night overflow wrap around modulo one week.

Edges (all directed, ``weight`` = seconds, attribute ``mode``):

* ``train``    -- a scheduled ride between consecutive stops of one trip.
                  Carries ``route``/``line``, ``trip``, ``direction``.
* ``wait``     -- staying at one platform until its next event (cyclic).
* ``transfer`` -- changing platforms within a station complex, or walking
                  between linked complexes (from ``transfers.txt``), honoring
                  ``min_transfer_time``. ``walk=True`` flags cross-complex hops.

Because nodes are platform-level, "visited a station" maps a platform to its
``parent_station`` (see :func:`station_of`).

Days are materialized from the GTFS ``calendar`` service patterns:
Weekday -> Mon..Fri, Saturday -> Sat, Sunday -> Sun.
"""
from __future__ import annotations

import argparse
import bisect
import pickle
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import pandas as pd

DAY = 86_400
WEEK = 7 * DAY

# Mon=0 .. Sun=6. Which weekday offsets each GTFS service pattern is active on.
SERVICE_TO_DAYS = {
    "Weekday": (0, 1, 2, 3, 4),
    "Saturday": (5,),
    "Sunday": (6,),
}

DEFAULT_GTFS = Path("data/gtfs")
DEFAULT_OUT = Path("data/graph.pkl")

# Routes bundled in gtfs_subway.zip but excluded by default. The Staten Island
# Railway (route_id "SI") has no track or in-system transfer link to the subway
# (ferry-only), so it is a disconnected island and is dropped for "subway only".
# NB: the Franklin Av Shuttle ("FS") also uses S0x stop ids but IS connected --
# always filter SIR by route_id, never by stop-id prefix.
DEFAULT_EXCLUDE_ROUTES = ("SI",)

# Reversal time at a terminal (arrival platform -> departure platform). The
# train reverses; the schedule gap to the next departure carries the layover, so
# a small nominal value suffices.
TERMINAL_REVERSAL_S = 60


def parse_gtfs_time(s: str) -> int:
    """``'HH:MM:SS'`` (HH may exceed 24) -> seconds since midnight."""
    h, m, sec = s.split(":")
    return int(h) * 3600 + int(m) * 60 + int(sec)


@dataclass
class Feed:
    stops: pd.DataFrame
    stop_times: pd.DataFrame
    trips: pd.DataFrame
    routes: pd.DataFrame
    transfers: pd.DataFrame

    @classmethod
    def load(cls, gtfs_dir: Path | str = DEFAULT_GTFS) -> "Feed":
        d = Path(gtfs_dir)
        return cls(
            stops=pd.read_csv(d / "stops.txt", dtype=str),
            stop_times=pd.read_csv(d / "stop_times.txt", dtype=str),
            trips=pd.read_csv(d / "trips.txt", dtype=str),
            routes=pd.read_csv(d / "routes.txt", dtype=str),
            transfers=pd.read_csv(d / "transfers.txt", dtype=str),
        )


def station_of(stop_id: str) -> str:
    """Map a platform stop id to its parent station id (drop trailing N/S)."""
    if stop_id and stop_id[-1] in "NS":
        return stop_id[:-1]
    return stop_id


class TimeExpandedGraphBuilder:
    def __init__(self, feed: Feed, days: tuple[int, ...] = tuple(range(7)),
                 routes: list[str] | None = None, add_transfers: bool = True,
                 exclude_routes: tuple[str, ...] = DEFAULT_EXCLUDE_ROUTES,
                 walk_links: list[tuple[str, str, int]] | None = None):
        self.feed = feed
        self.days = tuple(sorted(days))
        self.routes = set(routes) if routes else None
        self.exclude_routes = set(exclude_routes)
        self.add_transfers = add_transfers
        self.walk_links = walk_links or []

        self.G = nx.DiGraph()
        self._node_id: dict[tuple[str, int], int] = {}
        # platform stop_id -> sorted list of distinct event times present in graph
        self._timeline: dict[str, list[int]] = defaultdict(list)
        # parent station id -> set of platform stop_ids
        self._platforms: dict[str, set[str]] = defaultdict(set)

        self._route_line = feed.routes.set_index("route_id")["route_short_name"].to_dict()

    # -- node helpers -----------------------------------------------------
    def _node(self, stop: str, t: int) -> int:
        t %= WEEK
        key = (stop, t)
        nid = self._node_id.get(key)
        if nid is None:
            nid = len(self._node_id)
            self._node_id[key] = nid
            self.G.add_node(nid, stop=sys.intern(stop), t=t,
                            station=sys.intern(station_of(stop)))
            self._timeline[stop].append(t)
        return nid

    # -- build steps ------------------------------------------------------
    def build(self) -> nx.DiGraph:
        self._add_train_and_event_nodes()
        self._finalize_timelines()
        self._add_wait_edges()
        if self.add_transfers:
            self._add_transfer_edges()
        return self.G

    def _segments(self) -> pd.DataFrame:
        """One row per consecutive (from_stop -> to_stop) hop within a trip,
        with service pattern and route metadata attached."""
        st = self.feed.stop_times.copy()
        st["stop_sequence"] = st["stop_sequence"].astype(int)
        st["dep"] = st["departure_time"].map(parse_gtfs_time)
        st["arr"] = st["arrival_time"].map(parse_gtfs_time)

        trips = self.feed.trips
        if self.exclude_routes:
            trips = trips[~trips["route_id"].isin(self.exclude_routes)]
        if self.routes is not None:
            trips = trips[trips["route_id"].isin(self.routes)]
        st = st.merge(trips[["trip_id", "route_id", "service_id", "direction_id"]],
                      on="trip_id", how="inner")
        st = st[st["service_id"].isin(SERVICE_TO_DAYS)]
        st = st.sort_values(["trip_id", "stop_sequence"])

        nxt = st.groupby("trip_id", sort=False).shift(-1)
        seg = pd.DataFrame({
            "from_stop": st["stop_id"], "from_dep": st["dep"],
            "to_stop": nxt["stop_id"], "to_arr": nxt["arr"],
            "route": st["route_id"], "trip": st["trip_id"],
            "direction": st["direction_id"], "service": st["service_id"],
        })
        return seg.dropna(subset=["to_stop"])  # last stop of each trip has no next

    def _add_train_and_event_nodes(self) -> None:
        seg = self._segments()
        added = 0
        for r in seg.itertuples(index=False):
            line = self._route_line.get(r.route, r.route)
            for off in SERVICE_TO_DAYS[r.service]:
                if off not in self.days:
                    continue
                base = off * DAY
                u = self._node(r.from_stop, base + r.from_dep)
                v = self._node(r.to_stop, base + r.to_arr)
                w = (self.G.nodes[v]["t"] - self.G.nodes[u]["t"]) % WEEK
                self.G.add_edge(u, v, mode="train", weight=w, line=line,
                                route=sys.intern(r.route), trip=r.trip,
                                direction=r.direction)
                added += 1
        print(f"  train edges: {added}")

    def _finalize_timelines(self) -> None:
        for stop, times in self._timeline.items():
            uniq = sorted(set(times))
            self._timeline[stop] = uniq
            self._platforms[station_of(stop)].add(stop)

    def _add_wait_edges(self) -> None:
        added = 0
        for stop, times in self._timeline.items():
            n = len(times)
            if n < 2:
                continue
            for i in range(n):
                t0 = times[i]
                t1 = times[(i + 1) % n]
                w = (t1 - t0) % WEEK
                if w == 0:
                    continue
                self.G.add_edge(self._node(stop, t0), self._node(stop, t1),
                                mode="wait", weight=w)
                added += 1
        print(f"  wait edges: {added}")

    def _transfer_links(self) -> list[tuple[str, str, int, bool]]:
        """Unified (from_parent, to_parent, min_seconds, is_walk) links from
        GTFS transfers.txt, terminal reversals, and any out-of-system walk links.
        A walk link is just a transfer whose minimum time is the street-walk time."""
        links = []
        self_transfer = set()
        for tr in self.feed.transfers.itertuples(index=False):
            mtt = int(float(tr.min_transfer_time)) if pd.notna(tr.min_transfer_time) else 0
            links.append((tr.from_stop_id, tr.to_stop_id, mtt,
                          tr.from_stop_id != tr.to_stop_id))
            if tr.from_stop_id == tr.to_stop_id:
                self_transfer.add(tr.from_stop_id)

        # Terminal reversals: at a terminal the train reverses, so the arrival
        # platform connects to the departure platform. MTA's transfers.txt omits
        # this self-transfer for some terminals, leaving the arrival platform a
        # dead-end (e.g. South Ferry 142S). Add a reversal for any 2-platform
        # station that lacks a documented self-transfer.
        for parent, plats in self._platforms.items():
            if len(plats) >= 2 and parent not in self_transfer:
                links.append((parent, parent, TERMINAL_REVERSAL_S, False))

        for a, b, secs in self.walk_links:           # precise street walks
            links.append((a, b, int(secs), True))
        return links

    def _add_transfer_edges(self) -> None:
        """Connect each event at a source platform to the earliest feasible event
        at each target platform (>= the link's minimum time)."""
        added = 0
        for from_parent, to_parent, mtt, walk in self._transfer_links():
            src_plats = self._platforms.get(from_parent, set())
            dst_plats = self._platforms.get(to_parent, set())
            for sp in src_plats:
                src_times = self._timeline.get(sp, [])
                for dp in dst_plats:
                    if dp == sp:
                        continue
                    dst_times = self._timeline.get(dp, [])
                    if not dst_times:
                        continue
                    for t in src_times:
                        tt = self._first_after(dst_times, (t + mtt) % WEEK)
                        w = (tt - t) % WEEK
                        if w == 0:
                            continue
                        self.G.add_edge(self._node(sp, t), self._node(dp, tt),
                                        mode="transfer", weight=w, walk=walk,
                                        min_transfer_time=mtt)
                        added += 1
        print(f"  transfer edges: {added}")

    @staticmethod
    def _first_after(times: list[int], target: int) -> int:
        """Earliest time >= target in the cyclic list (wraps to times[0])."""
        i = bisect.bisect_left(times, target)
        return times[i] if i < len(times) else times[0]


def load_walk_links(path: Path | str) -> list[tuple[str, str, int]]:
    """Load (from_parent, to_parent, seconds) rows from a walk_transfers CSV."""
    df = pd.read_csv(path, dtype={"from_parent": str, "to_parent": str})
    return list(zip(df["from_parent"], df["to_parent"], df["seconds"].astype(int)))


def build_time_expanded_graph(gtfs_dir=DEFAULT_GTFS, days=tuple(range(7)),
                              routes=None, add_transfers=True,
                              exclude_routes=DEFAULT_EXCLUDE_ROUTES,
                              walk_links=None) -> nx.DiGraph:
    feed = Feed.load(gtfs_dir)
    G = TimeExpandedGraphBuilder(feed, days=days, routes=routes,
                                 add_transfers=add_transfers,
                                 exclude_routes=exclude_routes,
                                 walk_links=walk_links).build()
    print(f"  total: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Build the time-expanded subway graph.")
    p.add_argument("--gtfs", default=str(DEFAULT_GTFS))
    p.add_argument("--out", default=str(DEFAULT_OUT))
    p.add_argument("--days", default="0-6",
                   help="Day offsets to materialize, e.g. '0-6' (Mon-Sun) or '0'.")
    p.add_argument("--routes", default=None,
                   help="Comma-separated route_ids to restrict to (for testing).")
    p.add_argument("--exclude-routes", default=",".join(DEFAULT_EXCLUDE_ROUTES),
                   help="Comma-separated route_ids to drop (default: SI = Staten "
                        "Island Railway, which is disconnected from the subway).")
    p.add_argument("--no-transfers", action="store_true")
    p.add_argument("--walk-transfers", default=None,
                   help="CSV of precise walk links (from walk_transfers.py) to add.")
    args = p.parse_args(argv)

    if "-" in args.days:
        a, b = args.days.split("-")
        days = tuple(range(int(a), int(b) + 1))
    else:
        days = tuple(int(x) for x in args.days.split(","))
    routes = args.routes.split(",") if args.routes else None
    exclude = tuple(r for r in args.exclude_routes.split(",") if r)
    walk_links = load_walk_links(args.walk_transfers) if args.walk_transfers else None

    print(f"Building graph: days={days} routes={routes or 'ALL'} "
          f"exclude={exclude or 'none'} transfers={not args.no_transfers} "
          f"walk_links={len(walk_links) if walk_links else 0}")
    G = build_time_expanded_graph(args.gtfs, days=days, routes=routes,
                                  add_transfers=not args.no_transfers,
                                  exclude_routes=exclude, walk_links=walk_links)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as f:
        pickle.dump(G, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved -> {out} ({out.stat().st_size/1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

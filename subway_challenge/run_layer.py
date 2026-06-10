"""On-demand run-transfer layer for the time-expanded subway graph.

Instead of baking out-of-system run links into the graph as edges (which explodes
to ~93M edges at a 5 km radius), this layer computes run options *live* from a
solver's current node. A run can start at any time, so from node ``(stop, t)``
you may run to any nearby complex, arrive at ``t + run_seconds``, and board the
first train there -- no static edges, full radius coverage.

Usage in a solver::

    from subway_challenge.run_layer import RunLayer
    runs = RunLayer.from_graph(G, radius_m=5000)         # uses cached OSRM matrix
    for v, weight, info in runs.neighbors(node):         # graph edges + run options
        ...

``neighbors`` yields the node's ordinary graph successors *and* its run
successors uniformly, so a Dijkstra/A* can treat them identically. Use
``run_successors`` alone for just the runs.
"""
from __future__ import annotations

import bisect
from collections import defaultdict

from .build_graph import WEEK, station_of
from .walk_transfers import (DEFAULT_ACCESS_PENALTY_S, DEFAULT_PACE_MPS,
                             DEFAULT_RUN_RADIUS_M, complex_run_adjacency)


class RunLayer:
    def __init__(self, G, adjacency, parent_complex, complex_platforms, platform_index):
        self.G = G
        self.adjacency = adjacency                  # cid -> [(cid2, seconds, meters)]
        self.parent_complex = parent_complex        # parent stop id -> complex id
        self.complex_platforms = complex_platforms  # cid -> [platform stop ids in graph]
        self.platform_index = platform_index        # platform -> (times[], node_ids[])

    @classmethod
    def from_graph(cls, G, radius_m=DEFAULT_RUN_RADIUS_M, pace_mps=DEFAULT_PACE_MPS,
                   access_penalty_s=DEFAULT_ACCESS_PENALTY_S, max_seconds=None, **kw):
        adjacency, complexes = complex_run_adjacency(
            radius_m=radius_m, pace_mps=pace_mps, access_penalty_s=access_penalty_s,
            max_seconds=max_seconds, **kw)

        # Per-platform sorted (time, node_id) index, built once from the graph.
        tmp = defaultdict(list)
        for nid, d in G.nodes(data=True):
            tmp[d["stop"]].append((d["t"], nid))
        platform_index = {}
        for stop, lst in tmp.items():
            lst.sort()
            platform_index[stop] = ([t for t, _ in lst], [i for _, i in lst])

        parent_complex = {p: c.cid for c in complexes.values() for p in c.parents}
        complex_platforms = {}
        for cid, c in complexes.items():
            plats = [pl for parent in c.parents for pl in (parent + "N", parent + "S")
                     if pl in platform_index]
            complex_platforms[cid] = plats
        return cls(G, adjacency, parent_complex, complex_platforms, platform_index)

    def run_successors(self, node):
        """Yield ``(node_id, weight, info)`` run options from ``node``. ``weight``
        is total elapsed seconds (run + wait for the boarded train); ``info`` has
        mode, run_seconds, distance, and target complex."""
        d = self.G.nodes[node]
        t = d["t"]
        src = self.parent_complex.get(station_of(d["stop"]))
        if src is None:
            return
        for tgt, run_s, meters in self.adjacency.get(src, ()):
            arr = (t + run_s) % WEEK
            for plat in self.complex_platforms.get(tgt, ()):
                times, ids = self.platform_index[plat]
                i = bisect.bisect_left(times, arr)
                if i == len(times):              # wrap to first event next week
                    i = 0
                nid, tt = ids[i], times[i]
                yield nid, (tt - t) % WEEK, {
                    "mode": "run", "run_seconds": run_s, "meters": meters,
                    "to_complex": tgt, "board_platform": plat,
                }

    def neighbors(self, node):
        """Graph successors followed by run successors, as ``(node, weight, info)``."""
        for v in self.G.successors(node):
            e = self.G[node][v]
            yield v, e["weight"], e
        yield from self.run_successors(node)

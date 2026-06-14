# Subway Challenge — Codex Mission Brief

Goal: find and validate a route that visits all **472 official NYC subway
stations** in less than **22:14:10** (`80050` seconds), beating Kate Jones's
current world record. The current repo best is **24:45:00** (`89100` seconds),
so the working target is to save at least **02:30:51** and produce a validator
`RESULT` line with `elapsed_s < 80050`.

## Prime Directive

Do not merely improve the current best. The success condition is:

```text
RESULT valid=true stations=472/472 elapsed_s=<80050 elapsed=<22:14:10
```

Every serious candidate must be validated against the full solver:

```bash
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

`solutions/best.json` is the source of truth for the best known route. A faster
valid candidate should be recorded with `--record`; a route that looks good in a
planning graph but fails validation does not count.

## What Counts

* Visit identity is official MTA `Station ID`, not GTFS parent or complex.
* All **472** canonical stations must be covered.
* Legal moves are only those accepted by `solver.py`: train, wait, transfer, and
  run-layer transitions.
* Elapsed time is the sum of transition weights. Never score by raw timestamp
  subtraction without cyclic-week modulo handling.
* Start station, end station, and start time are free unless a specific
  experiment fixes them.

## Current Baseline

Known valid baseline:

```bash
venv/bin/python -m subway_challenge.solver best
# RESULT valid=true stations=472/472 elapsed_s=89100 elapsed=24:45:00 steps=664 modes={'train': 536, 'transfer': 43, 'wait': 56, 'run': 28}
```

Anatomy of the baseline:

* train time: about `16:05:30`
* run time charged: about `04:38:30`
* actual run distance: about `34.8 km`
* wait/transfer/run-wait overhead: nearly `5h`

This is not a "trim a few transfers" problem. To beat the record, the search
must find a different macro route basin with better branch ordering and better
schedule phase alignment.

## Strategy Learned So Far

The first-visit station-anchor LNS can find valid routes but gets trapped around
the `24:45:00` basin. Simple rotations, reversals, and local deletion of
zero-new rides do not produce a breakthrough.

Use these planning layers:

| Module | Purpose |
|---|---|
| `compact_graph.py` | Collapse passthrough corridors into decision-station macro-edges. |
| `compact_schedule.py` | Attach real GTFS scheduled corridor events to the compact graph. |
| `route_anatomy.py` | Diagnose runs, waits, repeated corridors, backtracks, and low-yield rides. |
| `route_surgery.py` | Try prefix-preserving local route surgeries and validate them. |
| `macro_search.py` | Search over scheduled compact corridor events, then realize on the full graph. |
| `search.py` | Existing exact LNS/seed sweep tools; still useful for polishing candidates. |

Key compact counts:

```text
compact graph: 212 decision stations, 260 passthrough stations collapsed
compact schedule: 734035 scheduled macro segments, 504 directed serviced edges
Pareto-dominated compact events: only about 0.5%
```

The low Pareto-dominance rate means naive "drop dominated lines/days" pruning is
not enough. Timing context matters.

## Next Best Approach

Build around **mandatory branch templates**, not unconstrained greedy coverage.
Every remote branch cluster must be covered intentionally and early enough that
the route does not leave a catastrophic repair tail.

High-value branch/template families include:

* Rockaways: Far Rockaway, Rockaway Park, Broad Channel, Lefferts.
* Canarsie/New Lots/Livonia area.
* Bronx IRT: Wakefield, Dyre, Pelham, Woodlawn, 148 St.
* Northern Manhattan/Bronx IND/BMT: Inwood, 207, Norwood, Bedford Park.
* Flushing/Astoria/Queens Boulevard/7.
* Coney Island and south Brooklyn branches.
* Myrtle/M/J/Z eastern branches.

A promising solver shape:

1. Enumerate branch templates with orientation choices.
2. Schedule each template using compact scheduled corridor events.
3. Connect templates with exact or bounded Dijkstra over the full graph/run
   layer.
4. Search over template order using beam/LNS/SA.
5. Realize the final station/platform path on `data/graph.pkl`.
6. Validate with `solver.py --record`.

## Search Discipline

When running a solving turn:

1. Start by validating the current best so the baseline is visible.
2. Generate one or more candidate routes in `solutions/`.
3. Validate every candidate with `--record`.
4. Surface every important `RESULT` line in the final response.
5. If no improvement is found, report the best failed/tied candidate and the
   concrete lesson learned.

Useful commands:

```bash
venv/bin/python -m subway_challenge.solver best
venv/bin/python -m subway_challenge.route_anatomy solutions/best.json --limit 12
venv/bin/python -m subway_challenge.compact_graph --out data/compact_graph.json
venv/bin/python -m subway_challenge.compact_schedule --out data/compact_schedule.json
venv/bin/python -m subway_challenge.search lns --seed-from solutions/best.json --terminal-runs --run-radius 2500
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

## Guardrails

* Keep `solver.py` as the validator/scorer contract.
* Do not trust compact or macro scores unless the realized route validates.
* Do not overfit to route count or number of runs; only elapsed time wins, but
  excessive running/switching is a smell when it creates repair tails.
* Keep generated probes in `solutions/` with descriptive names.
* Preserve `solutions/best.json` unless a candidate is valid and faster.

## Target

The mission is complete only when the repo contains a route with:

```text
valid=true
stations=472/472
elapsed_s < 80050
elapsed < 22:14:10
```

Until then, keep searching for a different macro basin. The problem is not
solved by polishing `24:45:00`; the world-record route requires a structural
breakthrough.

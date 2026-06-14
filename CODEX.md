# Subway Challenge — Codex Mission Brief

Goal: find a route that visits all **472 official NYC subway stations** in less
than **22:14:10** (`80050` seconds), beating Kate Jones's current world record.

The current repo best is **24:45:00** (`89100` seconds), so the required
improvement is at least **02:30:51**. Do not redefine success around a smaller
improvement. The mission is complete only when a route validates with:

```text
RESULT valid=true stations=472/472 elapsed_s=<80050
```

## Clean Start

Treat this repo as the starting point. Ignore previous experimental plans,
helper modules, and generated probe routes unless they are present in the clean
tree and still prove useful after inspection. Start from the code and data that
actually exist now.

The authoritative pieces are:

| Path | Role |
|---|---|
| `subway_challenge/solver.py` | Validator, scorer, and solution-format contract. |
| `subway_challenge/search.py` | Current LNS optimizer and best starting solver. |
| `subway_challenge/build_graph.py` | Builds the time-expanded GTFS graph. |
| `subway_challenge/run_layer.py` | Legal out-of-system run transitions. |
| `subway_challenge/walk_transfers.py` | Run-distance model and OSRM cache logic. |
| `subway_challenge/stations.py` | Official 472-station identity mapping. |
| `solutions/best.json` | Current best valid route. |

Do not assume any untracked or previously discussed helper exists. If a new
planning layer is needed, build it deliberately, keep it scoped, and validate
its output through `solver.py`.

## Non-Negotiable Validation

Every serious candidate must be checked with the full validator:

```bash
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

Always surface the `RESULT` line. A planning score, compact-graph score, visual
inspection, or plausible itinerary does not count unless the full validator says
the route is valid.

Current baseline command:

```bash
venv/bin/python -m subway_challenge.solver best
# expected current baseline: elapsed_s=89100 elapsed=24:45:00
```

`solutions/best.json` should change only through `--record` when a candidate is
valid and faster.

## Problem Rules

* Visit identity is official MTA `Station ID`; all **472** canonical stations
  must be covered.
* Legal moves are train, wait, transfer, and run transitions accepted by
  `solver.py`.
* Score is total elapsed seconds: sum transition weights, not raw timestamp
  differences.
* The timetable is cyclic over a week. Any custom elapsed calculation must
  handle week wrap correctly.
* Start station, end station, and start time are free unless an experiment
  explicitly fixes them.

## Search Principles

The current best is already far from a naive route. To beat `22:14:10`, expect a
structural route change, not just a small cleanup.

Pursue ideas that can plausibly save hours:

* Start-time and day-type sweeps.
* Alternative starts/finishes, especially human-plausible terminal choices.
* Branch-order search instead of only first-visit station-order tweaks.
* Explicit handling of remote branch clusters so they are not left for a late
  repair tail.
* Run-aware but not run-addicted routing; excessive running can mask bad macro
  order.
* Schedule-aware construction and validation on the real GTFS graph.

Avoid spending too much time on ideas that can only save seconds:

* Pure formatting/refactoring.
* Local route polishing without evidence it changes the macro route.
* Dropping lines/days globally without proving dominance in the schedule
  context.

## Working Loop

For each solving turn:

1. Inspect current state rather than relying on memory.
2. Validate the current best.
3. Make a concrete solver/search improvement or produce a candidate.
4. Validate candidates with `--record`.
5. Report the best `RESULT` line and the lesson learned.

Useful commands:

```bash
venv/bin/python -m subway_challenge.solver best
venv/bin/python -m subway_challenge.search lns --seed-from solutions/best.json --terminal-runs --run-radius 2500
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

## Completion Target

The goal remains open until the repo contains a route proven by `solver.py` with:

```text
valid=true
stations=472/472
elapsed_s < 80050
elapsed < 22:14:10
```

Anything else is progress, not completion.

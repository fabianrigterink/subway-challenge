# Subway Challenge - Codex Charter

This is the operating brief for Codex in this repository. Treat the checkout as
a clean starting point for solving the NYC Subway Challenge with the MTA
schedule data already present in the repo.

## Mission

Find and validate a route that visits all **472 official NYC subway stations**
in strictly less than **22:14:10**.

That is the bar. Improving the repository baseline is useful, but the mission is
not complete unless the validated elapsed time is below Kate Jones's `22:14:10`
benchmark:

```text
valid=true
stations=472/472
elapsed_s < 80050
elapsed < 22:14:10
```

A tie does not count. A heuristic estimate does not count. A partial route does
not count. Only the full validator decides.

## Clean-Start Rule

Start from scratch on every serious solving pass.

Ignore prior conversations, remembered experiments, stale notes, untracked probe
files, and old helper scripts unless you deliberately inspect and revalidate
them. The repository and validator are the source of truth, not the search
history.

Begin by orienting from the current checkout:

```bash
git status --short
git ls-files
venv/bin/python -m subway_challenge.solver best
```

If the working tree contains dirty or untracked files from earlier attempts,
leave them alone unless they are directly relevant. Treat them as scratch until
they have been read, understood, and validated.

## Source Of Truth

Read these files before changing solver behavior:

| Path | Purpose |
|---|---|
| `README.md` | Project overview, baseline, data pipeline, and known commands. |
| `CLAUDE.md` | Existing solver contract and validation expectations. |
| `subway_challenge/solver.py` | Validator, scorer, legal moves, and solution format. |
| `subway_challenge/search.py` | Main search and route-construction entry point. |
| `subway_challenge/build_graph.py` | GTFS-derived time-expanded graph construction. |
| `subway_challenge/run_layer.py` | Legal out-of-system run transitions. |
| `subway_challenge/walk_transfers.py` | Run-distance and OSRM model. |
| `subway_challenge/stations.py` | Official 472-station identity mapping. |
| `solutions/best.json` | Current best route, trusted only after validation. |

Do not weaken the validator, station identity mapping, legal move model, or
elapsed-time scoring to make a route pass. Improve the route, not the rules.

## Validation Contract

Every meaningful candidate must be validated on the full graph:

```bash
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

Always surface the complete `RESULT` line. `--record` may update
`solutions/best.json`, but only when the candidate is valid and faster than the
current recorded best.

Compact graphs, lower bounds, local scores, visual inspection, and plausible
itineraries are allowed as diagnostics. They are not evidence of success until
the full validator reports a valid 472-station route.

## Problem Rules

* Coverage is by official MTA `Station ID`; all **472** canonical stations must
  be visited.
* Legal transitions are exactly the moves accepted by `solver.py`: `train`,
  `wait`, `transfer`, and `run`.
* Cost is total elapsed seconds, computed from transition weights.
* The timetable is cyclic over a week; never score by subtracting raw timestamps
  without handling wraparound.
* Start station, end station, service day, and start time are free unless an
  experiment explicitly fixes them.
* Use the actual GTFS-derived schedules in this repo. Simplified or synthetic
  models can guide search, but final candidates must realize on the full graph.

## Search Priorities

The current repository best is around **24:24:30**, while the target is under
**22:14:10**.
That gap is too large for cosmetic polishing. Prefer structural search moves
that can change the route shape.

Promising directions:

* Simplify the network to true decision points while preserving recoverable
  platform-level paths through pass-through stations.
* Compare services, lines, branches, time windows, and day types for dominance
  before pruning them.
* Penalize or temporarily forbid long waits, excessive transfers, and excessive
  out-of-system running during construction, then validate final routes under
  the unmodified elapsed-time objective.
* Search over start station, finish station, start time, weekday/weekend phase,
  branch order, and trunk traversal order.
* Build compact or macro models for idea generation, then realize candidates on
  the full time-expanded graph.
* Track lower bounds, bottlenecks, and failed assumptions so discarded
  directions still teach something.

Avoid spending a whole session only on small local edits around the incumbent
unless the experiment is designed to test a clear hypothesis.

## Working Loop

For each solving session:

1. Inspect the repo state.
2. Validate the current best and record its `RESULT`.
3. Choose one concrete modeling, pruning, construction, or improvement-search
   hypothesis.
4. Implement the smallest useful change or run the smallest decisive experiment.
5. Write candidate routes under `solutions/`.
6. Validate candidates with `--record`.
7. Report the best `RESULT` line and the lesson learned.

Useful commands:

```bash
venv/bin/python -m subway_challenge.solver best
venv/bin/python -m subway_challenge.search lns --seed-from solutions/best.json --terminal-runs --run-radius 2500
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

## Completion

Keep going until the repository contains a candidate proven by the full validator
to satisfy:

```text
valid=true
stations=472/472
elapsed_s < 80050
elapsed < 22:14:10
```

Anything else is a baseline improvement, an experiment, or a useful clue. It is
not the finish line.

# solutions/

Candidate and best-known routes for the Subway Challenge.

* `best.json` — the best valid route so far (written by `solver.py validate
  --record` when a faster valid route is found). Tracked in git.
* `candidate.json` (or any name) — work-in-progress solutions the solver writes
  and validates each turn.

Format and rules: see [`../CLAUDE.md`](../CLAUDE.md) and
[`../subway_challenge/solver.py`](../subway_challenge/solver.py).

```json
{"meta": {"radius_m": 5000, "notes": "greedy seed"},
 "path": [["127S", 28860], ["125S", 28950]]}
```

Validate / score:

```bash
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
venv/bin/python -m subway_challenge.solver best
```

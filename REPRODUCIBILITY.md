# Reproducibility

This repository does not commit the generated `data/` cache. It is about
886 MB locally and includes the GTFS-derived time-expanded graph, compact
analysis artifacts, and the OSRM street-distance cache.

## Exact Local Validation

The committed `solutions/best.json` was validated against this local data cache:

| Artifact | SHA-256 |
|---|---|
| `data/graph.pkl` | `006b22ef794d366ec6766b0481191d397765f3d86fce6dff5cc8fa3886025038` |
| `data/walk/osrm_matrix.json` | `44ccb659ea5d2a0915ebc7f707bc42153aa209315b0a072d0cc55d8d24056320` |
| `data/gtfs/feed_info.txt` | `e1bbac08ffdc683a5d11d9aa332885681a79e4215a18d119ac6bb8dc0346bbe0` |

GTFS feed metadata:

```text
feed_publisher_name: MTA New York City Transit
feed_start_date: 20260526
feed_end_date: 20260907
feed_version: 20260526
```

Check the current local cache:

```bash
shasum -a 256 data/graph.pkl data/walk/osrm_matrix.json data/gtfs/feed_info.txt
venv/bin/python -m subway_challenge.solver best
```

Expected route result:

```text
RESULT valid=true stations=472/472 elapsed_s=87870 elapsed=24:24:30
```

## Fresh Rebuild

A fresh clone can rebuild the ignored data artifacts:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python -m subway_challenge.download_gtfs
python -m subway_challenge.build_graph
python -m subway_challenge.solver best
```

Because `download_gtfs` fetches the current public MTA feed, a future rebuild may
use different schedules. In that case, the code remains reproducible, but exact
route replay is only guaranteed if the rebuilt graph matches the hashes above.

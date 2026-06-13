"""Render a Subway Challenge solution as an interactive HTML map replay.

Reads a solution file (``solutions/*.json``) plus the already-built data layer
(the time-expanded graph + GTFS for track geometry and official MTA line colors)
and emits ONE ``.html`` file: a zoomable Leaflet map (Carto basemap) with the
subway network in real MTA colors, the solution route that solidifies as it is
ridden, a gliding marker, a synced itinerary, a live time breakdown, and a
stations-visited counter that lights up dots in their line color as covered.

Needs a network connection when viewed (Leaflet + map tiles load from a CDN).

Usage::

    python -m subway_challenge.visualize_map                        # -> best_replay_map.html
    python -m subway_challenge.visualize_map solutions/best.json -o my_map.html
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

from .build_graph import WEEK, station_of
from .solver import _load_solution, hms, load_problem

GTFS = Path("data/gtfs")
DAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]


def _load_stops() -> dict:
    """stop_id -> {lat, lon, name, parent}."""
    stops = {}
    with open(GTFS / "stops.txt", newline="") as f:
        for row in csv.DictReader(f):
            try:
                lat, lon = float(row["stop_lat"]), float(row["stop_lon"])
            except (ValueError, KeyError):
                continue
            stops[row["stop_id"]] = {
                "lat": lat, "lon": lon,
                "name": row["stop_name"],
                "parent": row.get("parent_station") or row["stop_id"],
            }
    return stops


def _load_route_colors() -> dict:
    """route_id -> {color, text} as #rrggbb."""
    colors = {}
    with open(GTFS / "routes.txt", newline="") as f:
        for row in csv.DictReader(f):
            c = (row.get("route_color") or "").strip()
            t = (row.get("route_text_color") or "").strip()
            colors[row["route_id"]] = {
                "color": f"#{c}" if c else "#6d6e71",
                "text": f"#{t}" if t else "#ffffff",
            }
    return colors


def _perp(p, a, b):
    """Planar perpendicular distance of p from segment a-b (lon scaled by cos lat)."""
    kx = math.cos(math.radians(a[0]))
    px, py = (p[1] - a[1]) * kx, p[0] - a[0]
    bx, by = (b[1] - a[1]) * kx, b[0] - a[0]
    L = math.hypot(bx, by)
    return math.hypot(px, py) if L == 0 else abs(px * by - py * bx) / L


def _rdp(pts, eps):
    """Ramer-Douglas-Peucker polyline simplification."""
    if len(pts) < 3:
        return pts
    a, b, dmax, idx = pts[0], pts[-1], 0.0, 0
    for i in range(1, len(pts) - 1):
        d = _perp(pts[i], a, b)
        if d > dmax:
            dmax, idx = d, i
    if dmax > eps:
        return _rdp(pts[:idx + 1], eps)[:-1] + _rdp(pts[idx:], eps)
    return [a, b]


def _load_shapes(route_colors: dict, shape_pts: dict, eps: float = 6e-5, exclude=("SI",)) -> list:
    """Real track geometry (from ``shape_pts``) -> simplified polylines colored by
    route, with reversed-direction duplicates dropped. Gives smooth, curved lines
    that follow the actual track instead of straight station-to-station chords.
    Routes in ``exclude`` (Staten Island Railway by default) are dropped."""
    sys.setrecursionlimit(10000)
    sroute = {}
    with open(GTFS / "trips.txt", newline="") as f:
        for row in csv.DictReader(f):
            if row.get("shape_id"):
                sroute[row["shape_id"]] = row["route_id"]
    out, seen = [], set()
    for sid, seq in shape_pts.items():
        if sroute.get(sid) in exclude:        # skip Staten Island Railway
            continue
        poly = _rdp(seq, eps)
        if len(poly) < 2:
            continue
        key = tuple((round(la, 5), round(lo, 5)) for la, lo in poly)
        if key in seen or key[::-1] in seen:      # drop the opposite-direction twin
            continue
        seen.add(key)
        col = route_colors.get(sroute.get(sid, ""), {}).get("color", "#6d6e71")
        out.append({"pts": [[round(la, 5), round(lo, 5)] for la, lo in poly], "c": col})
    return out


def _load_shape_points() -> dict:
    """Raw ordered track geometry: shape_id -> [(lat, lon), ...]."""
    pts = {}
    with open(GTFS / "shapes.txt", newline="") as f:
        for row in csv.DictReader(f):
            pts.setdefault(row["shape_id"], []).append(
                (int(row["shape_pt_sequence"]), float(row["shape_pt_lat"]), float(row["shape_pt_lon"])))
    return {sid: [(la, lo) for _, la, lo in sorted(s)] for sid, s in pts.items()}


def _trip_shapes() -> dict:
    """trip_id -> shape_id."""
    m = {}
    with open(GTFS / "trips.txt", newline="") as f:
        for row in csv.DictReader(f):
            m[row["trip_id"]] = row.get("shape_id", "")
    return m


def _subpath(pts, a, b, kx, eps=3e-5):
    """The slice of shape ``pts`` between the vertices nearest a and b -- i.e. the
    REAL track ridden between two stops, so it lies on the backdrop line."""
    if not pts:
        return None

    def near(p):
        return min(range(len(pts)),
                   key=lambda i: (pts[i][0] - p[0]) ** 2 + ((pts[i][1] - p[1]) * kx) ** 2)

    ia, ib = near(a), near(b)
    if ia == ib:
        return None
    seg = pts[ia:ib + 1] if ia < ib else pts[ib:ia + 1][::-1]
    return [[round(la, 5), round(lo, 5)] for la, lo in _rdp(seg, eps)]


def _coords(stops: dict, stop_id: str):
    """(lat, lon) for a platform or parent stop, falling back to the parent."""
    s = stops.get(stop_id) or stops.get(station_of(stop_id))
    return (s["lat"], s["lon"]) if s else None


def build_payload(solution_file: Path, radius_m: int | None = None) -> dict:
    path, meta = _load_solution(solution_file)
    prob = load_problem(radius_m=meta.get("radius_m", radius_m or 5000))
    G = prob.G
    stops = _load_stops()
    route_colors = _load_route_colors()

    used_routes = set()

    # --- station dots (GTFS parents that appear in the network) ---
    parents = {G.nodes[x]["station"] for x in G.nodes}
    dots = []
    dot_index = {}
    for sid in sorted(parents):
        c = _coords(stops, sid)
        if not c:
            continue
        dot_index[sid] = len(dots)
        dots.append({
            "lat": c[0], "lon": c[1],
            "name": stops.get(sid, {}).get("name", sid),
            "off": prob.stations.resolve(sid),   # official Station ID (or None)
            "vis": None,                          # elapsed-s when first visited
        })

    # --- the taken path: waypoints with cumulative elapsed + per-segment mode ---
    shape_pts = _load_shape_points()
    trip2shape = _trip_shapes()
    kx_proj = math.cos(math.radians(40.7))
    nodes = [prob.node_of(stop, int(t)) for stop, t in path]
    waypoints = []
    elapsed = 0.0
    start_t = int(path[0][1])
    visited_official = {}          # official id -> first elapsed
    cur_color = None               # color of the line currently being ridden (visit tint)
    for i, nid in enumerate(nodes):
        nd = G.nodes[nid]
        c = _coords(stops, nd["stop"]) or (0.0, 0.0)
        seg_mode = seg_route = seg_color = None
        seg_run_s = seg_run_m = 0
        seg_track = None
        if i > 0:
            w, mode = prob.transition(nodes[i - 1], nid)
            elapsed += w or 0
            seg_mode = mode
            if mode == "train":
                seg_route = G[nodes[i - 1]][nid].get("route", "")
                seg_color = route_colors.get(seg_route, {}).get("color", "#6d6e71")
                used_routes.add(seg_route)
                cur_color = seg_color
                sid = trip2shape.get(G[nodes[i - 1]][nid].get("trip"))
                ca = _coords(stops, G.nodes[nodes[i - 1]]["stop"])
                cb = _coords(stops, nd["stop"])
                if sid in shape_pts and ca and cb:
                    seg_track = _subpath(shape_pts[sid], ca, cb, kx_proj)
            elif mode == "run":
                best = None                      # match transition()'s min-weight run
                for rv, rw, info in prob.runs.run_successors(nodes[i - 1]):
                    if rv == nid and (best is None or rw < best[0]):
                        best = (rw, info)
                if best:
                    seg_run_s = int(best[1]["run_seconds"])
                    seg_run_m = best[1]["meters"]
        waypoints.append({
            "lat": c[0], "lon": c[1],
            "stop": nd["stop"],
            "name": stops.get(nd["stop"], {}).get("name")
                    or stops.get(station_of(nd["stop"]), {}).get("name", nd["stop"]),
            "e": round(elapsed),
            "mode": seg_mode,        # mode of the segment ARRIVING at this waypoint
            "route": seg_route,
            "color": seg_color,
            "rs": seg_run_s,         # actual running seconds (run segs; rest is waiting)
            "rm": round(seg_run_m),  # run distance in meters
            "track": seg_track,      # real ridden track (train segs) -> matches backdrop
        })
        off = prob.stations.resolve(nd["stop"])
        if off and off not in visited_official:
            visited_official[off] = round(elapsed)
            idx = dot_index.get(nd["station"])
            if idx is not None:
                dots[idx]["vis"] = round(elapsed)
                dots[idx]["vc"] = cur_color        # line color it was visited on

    # backfill the start station (visited before any train) with the first line color
    first_color = next((w["color"] for w in waypoints if w["color"]), "#3b3b3b")
    for d in dots:
        if d.get("vis") is not None and not d.get("vc"):
            d["vc"] = first_color

    legend = sorted(
        ({"route": r, "color": route_colors.get(r, {}).get("color", "#6d6e71"),
          "text": route_colors.get(r, {}).get("text", "#fff")}
         for r in used_routes if r),
        key=lambda x: x["route"],
    )

    itinerary = _itinerary(waypoints)

    return {
        "meta": meta,
        "file": str(solution_file),
        "total_s": round(elapsed),
        "total_hms": hms(elapsed),
        "start_t": start_t,
        "n_official": len(prob.canonical),
        "n_visited": len(set(visited_official) & prob.canonical),
        "shapes": _load_shapes(route_colors, shape_pts),
        "dots": dots,
        "waypoints": waypoints,
        "itinerary": itinerary,
        "legend": legend,
        "days": DAYS,
        "week": WEEK,
    }


def _itinerary(wp: list) -> list:
    """Collapse waypoints into a human trip log: consecutive same-line hops ->
    one 'ride', plus 'run' and 'connect' (transfer/wait) entries. Each entry
    carries its elapsed window so the replay can sync/highlight it."""
    n = len(wp)
    smode = lambda s: wp[s + 1]["mode"]          # mode of segment s (arriving at wp[s+1])
    sroute = lambda s: wp[s + 1]["route"]
    out, s = [], 0
    while s < n - 1:
        m = smode(s)
        if m == "train":
            r, e = sroute(s), s
            while e + 1 < n - 1 and smode(e + 1) == "train" and sroute(e + 1) == r:
                e += 1
            out.append({"type": "ride", "line": r, "color": wp[e + 1]["color"] or "#111",
                        "frm": wp[s]["name"], "to": wp[e + 1]["name"], "stops": e - s + 1,
                        "e0": wp[s]["e"], "e1": wp[e + 1]["e"]})
            s = e + 1
        elif m == "run":
            out.append({"type": "run", "color": "#e02b2b",
                        "frm": wp[s]["name"], "to": wp[s + 1]["name"],
                        "e0": wp[s]["e"], "e1": wp[s + 1]["e"]})
            s += 1
        else:  # transfer / wait bundle
            e = s
            while e + 1 < n - 1 and smode(e + 1) in ("transfer", "wait"):
                e += 1
            has_tr = any(smode(k) == "transfer" for k in range(s, e + 1))
            out.append({"type": "connect", "color": "#888",
                        "label": "Transfer" if has_tr else "Wait",
                        "frm": wp[s]["name"], "to": wp[e + 1]["name"],
                        "e0": wp[s]["e"], "e1": wp[e + 1]["e"]})
            s = e + 1
    return out


# --- Leaflet variant: overlay the route on a real zoomable geographic basemap ---
_LEAFLET_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Subway Challenge — replay</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
  :root { --ink:#1a1a1a; --panel:#ffffff; }
  * { box-sizing:border-box; }
  html,body { margin:0; height:100%; font-family:'Helvetica Neue',Helvetica,Arial,sans-serif; color:var(--ink); }
  #app { display:flex; flex-direction:column; height:100%; }
  #top { flex:1; display:flex; min-height:0; }
  #map { flex:1; min-height:0; background:#e9e5dc; }
  #side { width:330px; flex:none; overflow-y:auto; background:#fff; border-left:1px solid #ddd; }
  #side h3 { margin:0; padding:11px 14px; font-size:13px; border-bottom:1px solid #eee; position:sticky; top:0; background:#fff; z-index:1; }
  #side h3 span { font-weight:400; color:#888; }
  .it { display:flex; gap:10px; padding:8px 12px; border-bottom:1px solid #f2f2f2; cursor:pointer; border-left:3px solid transparent; }
  .it:hover { background:#fafafa; }
  .it.cur { background:#fff7e6; border-left-color:#222; }
  .it .ic { flex:none; width:26px; height:26px; border-radius:50%; display:flex; align-items:center; justify-content:center; font-weight:700; font-size:12px; color:#fff; margin-top:1px; }
  .it .t { font-size:13px; font-weight:600; }
  .it .s { font-size:11px; color:#666; margin-top:1px; font-variant-numeric:tabular-nums; }
  @media (max-width:760px){ #side{ display:none; } }
  .float { position:absolute; z-index:1200; background:rgba(255,255,255,.93); border:1px solid #ddd; border-radius:8px; box-shadow:0 2px 6px rgba(0,0,0,.12); font-size:12px; }
  #title { top:10px; left:54px; padding:8px 11px; } #title b { font-size:14px; }
  #legend { top:10px; right:10px; padding:8px 10px; max-width:172px; }
  #breakdown { bottom:16px; left:12px; padding:10px 12px; width:220px; }
  #breakdown .bd-h { font-weight:700; margin-bottom:7px; }
  #breakdown .bd-h span { font-weight:400; color:#888; }
  .bd-row { display:flex; justify-content:space-between; font-size:12px; margin-top:7px; }
  .bd-l { display:flex; align-items:center; gap:6px; }
  .bd-c { width:10px; height:10px; border-radius:2px; }
  .bd-v { font-variant-numeric:tabular-nums; color:#333; }
  .bd-track { height:5px; background:#eee; border-radius:3px; margin-top:3px; overflow:hidden; }
  .bd-track i { display:block; height:100%; width:0%; transition:width .12s linear; }
  .bd-dist { margin-top:10px; padding-top:8px; border-top:1px solid #eee; font-size:12px; }
  .bd-dist b { font-variant-numeric:tabular-nums; font-size:14px; }
  #legend .chip { display:inline-flex; align-items:center; justify-content:center; width:18px; height:18px; border-radius:50%; font-weight:700; margin:2px; font-size:11px; }
  #bar { background:var(--panel); border-top:1px solid #ddd; padding:10px 16px; box-shadow:0 -2px 8px rgba(0,0,0,.08); }
  #readout { display:flex; gap:22px; align-items:baseline; flex-wrap:wrap; margin-bottom:8px; }
  #clock { font-size:22px; font-weight:700; font-variant-numeric:tabular-nums; }
  #day { font-size:13px; color:#666; font-weight:600; }
  .stat { font-size:13px; color:#444; } .stat b { font-size:16px; font-variant-numeric:tabular-nums; }
  #status { font-weight:600; }
  #ctrl { display:flex; gap:10px; align-items:center; }
  #slider { flex:1; height:6px; }
  button { font:inherit; border:1px solid #bbb; background:#fff; border-radius:7px; padding:5px 11px; cursor:pointer; }
  button:hover { background:#f0f0f0; } button.on { background:var(--ink); color:#fff; border-color:var(--ink); }
  .marker-halo { border-radius:50%; animation:pulse 1.5s ease-out infinite; }
  @keyframes pulse { 0%{box-shadow:0 0 0 0 var(--pc, rgba(0,0,0,.5))} 100%{box-shadow:0 0 0 16px rgba(0,0,0,0)} }
</style>
</head>
<body>
<div id="app">
  <div id="top">
    <div id="map">
      <div id="title" class="float"></div>
      <div id="legend" class="float"></div>
      <div id="breakdown" class="float">
        <div class="bd-h">Time so far <span id="bdTotal"></span></div>
        <div class="bd-row"><span class="bd-l"><span class="bd-c" style="background:#2f7ed8"></span>Riding</span><span class="bd-v" id="bdRide">—</span></div>
        <div class="bd-track"><i id="barRide" style="background:#2f7ed8"></i></div>
        <div class="bd-row"><span class="bd-l"><span class="bd-c" style="background:#e02b2b"></span>Running</span><span class="bd-v" id="bdRun">—</span></div>
        <div class="bd-track"><i id="barRun" style="background:#e02b2b"></i></div>
        <div class="bd-row"><span class="bd-l"><span class="bd-c" style="background:#f5a623"></span>Waiting</span><span class="bd-v" id="bdWait">—</span></div>
        <div class="bd-track"><i id="barWait" style="background:#f5a623"></i></div>
        <div class="bd-dist" id="bdDist"></div>
      </div>
    </div>
    <div id="side"><h3>Itinerary <span id="itcount"></span></h3><div id="itlist"></div></div>
  </div>
  <div id="bar">
    <div id="readout">
      <span><span id="clock">--:--:--</span> <span id="day"></span></span>
      <span class="stat">elapsed <b id="elapsed">00:00:00</b> / <span id="total"></span></span>
      <span class="stat">stations <b id="visited">0</b>/<span id="nstations"></span></span>
      <span class="stat" id="status">—</span>
    </div>
    <div id="ctrl">
      <button id="play">▶ Play</button>
      <input id="slider" type="range" min="0" max="1000" value="0" step="1">
      <span>speed</span>
      <button class="speed" data-mult="200">200×</button>
      <button class="speed on" data-mult="800">800×</button>
      <button class="speed" data-mult="3000">3000×</button>
    </div>
  </div>
</div>
<script>
const DATA = /*__DATA__*/;
const map = L.map('map', {zoomControl:true, minZoom:9, maxZoom:17, preferCanvas:false});
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png', {
  subdomains:'abcd', maxZoom:19, attribution:'&copy; <a href="https://openstreetmap.org/copyright">OpenStreetMap</a> &copy; <a href="https://carto.com/attributions">CARTO</a>'
}).addTo(map);
map.fitBounds(L.latLngBounds(DATA.dots.map(d=>[d.lat,d.lon])).pad(0.04));

// stacking order: skeleton < route < stations < moving marker
map.createPane('pSkel');    map.getPane('pSkel').style.zIndex    = 405;
map.createPane('pRoute');   map.getPane('pRoute').style.zIndex   = 412;
map.createPane('pStation'); map.getPane('pStation').style.zIndex = 420;
map.createPane('pTop');     map.getPane('pTop').style.zIndex     = 480;
const skelR  = L.canvas({pane:'pSkel',    padding:0.5});
const routeR = L.svg   ({pane:'pRoute',   padding:0.5});
const dotR   = L.canvas({pane:'pStation', padding:0.5});
const topR   = L.svg   ({pane:'pTop',     padding:0.5});

// Chaikin corner-cutting: rounds a polyline by replacing each corner with two
// points at 1/4 and 3/4 of its edges (endpoints fixed). More iterations = rounder.
function chaikin(pts, it){
  let p=pts;
  for(let n=0;n<it;n++){
    if(p.length<3) break;
    const q=[p[0]];
    for(let i=0;i<p.length-1;i++){ const a=p[i], b=p[i+1];
      q.push([a[0]*0.75+b[0]*0.25, a[1]*0.75+b[1]*0.25]);
      q.push([a[0]*0.25+b[0]*0.75, a[1]*0.25+b[1]*0.75]); }
    q.push(p[p.length-1]); p=q;
  }
  return p;
}

// pre-lighten a hex color toward white (used so the faint backdrop is drawn at
// FULL opacity -- overlapping same-color services then overdraw the identical
// pale tint instead of stacking translucency and darkening the trunk).
function lighten(hex,t){ const c=hex.replace('#',''); const r=parseInt(c.slice(0,2),16),g=parseInt(c.slice(2,4),16),b=parseInt(c.slice(4,6),16);
  const m=v=>Math.round(v+(255-v)*t); return `rgb(${m(r)},${m(g)},${m(b)})`; }
function rgba(hex,a){ const c=hex.replace('#',''); return `rgba(${parseInt(c.slice(0,2),16)},${parseInt(c.slice(2,4),16)},${parseInt(c.slice(4,6),16)},${a})`; }

// --- network from real track geometry (curved, MTA-colored, corner-rounded) ---
for (const s of DATA.shapes)
  L.polyline(chaikin(s.pts,3), {renderer:skelR, color:lighten(s.c,0.72), weight:2.6, opacity:1, interactive:false,
    lineJoin:'round', lineCap:'round'}).addTo(map);

// --- station dots (drawn above the route so visited ones stay visible) ---
// unvisited: small white tick w/ thin grey ring; visited: solid line-color bullet, no border
// stations carry the progress: not-yet-visited = small hollow grey tick; visited = bold line-color dot
const DOT_UNVIS = {stroke:true, color:'#c7cbd0', weight:1, fillColor:'#ffffff', fillOpacity:1, radius:2.2};
const DOT_VIS   = {stroke:true, color:'#ffffff', weight:1.5, fillOpacity:1, radius:4.4};   // fillColor (line) set per dot
const dotLayers = DATA.dots.map(d =>
  L.circleMarker([d.lat,d.lon], {renderer:dotR, ...DOT_UNVIS}).bindTooltip(d.name).addTo(map));

// --- solution route: pre-drawn per-segment, revealed as the marker passes ---
const WP = DATA.waypoints, LL = WP.map(w=>[w.lat,w.lon]);
function segColor(i){ const b=WP[i+1]; return b.mode==='train' ? (b.color||'#111') : (b.mode==='run' ? '#e02b2b' : '#888'); }

// --- time-breakdown prefix sums: riding / running / waiting (+ run distance) ---
// a 'run' segment = run_seconds of running, then waiting for the boarded train.
const pRide=new Float64Array(WP.length), pRun=new Float64Array(WP.length),
      pWait=new Float64Array(WP.length), pDist=new Float64Array(WP.length);
for (let s=0;s<WP.length-1;s++){
  const b=WP[s+1], w=b.e-WP[s].e; let dr=0,dn=0,dw=0,dd=0;
  if (b.mode==='train') dr=w;
  else if (b.mode==='run'){ dn=Math.min(b.rs,w); dw=Math.max(0,w-b.rs); dd=b.rm; }
  else dw=w;                                   // wait / transfer
  pRide[s+1]=pRide[s]+dr; pRun[s+1]=pRun[s]+dn; pWait[s+1]=pWait[s]+dw; pDist[s+1]=pDist[s]+dd;
}
function updateBreakdown(e){
  const i=findSeg(e), A=WP[i], B=WP[i+1], part=Math.max(0,e-A.e);
  let ride=pRide[i], run=pRun[i], wait=pWait[i], dist=pDist[i];
  if (B.mode==='train') ride+=part;
  else if (B.mode==='run'){ const rs=B.rs; run+=Math.min(part,rs); wait+=Math.max(0,part-rs);
    dist += B.rm * (rs>0 ? Math.min(part,rs)/rs : (part>0?1:0)); }
  else wait+=part;
  const tot=Math.max(1,ride+run+wait), pct=x=>Math.round(x/tot*100);
  bdTotal.textContent='· '+fmt(tot);
  bdRide.textContent=`${fmt(ride)} · ${pct(ride)}%`; barRide.style.width=(ride/tot*100)+'%';
  bdRun.textContent =`${fmt(run)} · ${pct(run)}%`;   barRun.style.width =(run/tot*100)+'%';
  bdWait.textContent=`${fmt(wait)} · ${pct(wait)}%`; barWait.style.width=(wait/tot*100)+'%';
  const km=dist/1000;
  bdDist.innerHTML=`🏃 Run distance <b>${km.toFixed(2)} km</b> <span style="color:#888">(${(km*0.621371).toFixed(2)} mi)</span>`;
}
const bdTotal=document.getElementById('bdTotal'), bdRide=document.getElementById('bdRide'),
      bdRun=document.getElementById('bdRun'), bdWait=document.getElementById('bdWait'),
      bdDist=document.getElementById('bdDist'), barRide=document.getElementById('barRide'),
      barRun=document.getElementById('barRun'), barWait=document.getElementById('barWait');

// --- ridden geometry follows the REAL track (same shapes/rounding as the backdrop) ---
const segGeom = WP.slice(0,-1).map((_,i)=>{
  const B=WP[i+1];
  return (B.mode==='train' && B.track && B.track.length>1) ? chaikin(B.track,3) : [LL[i],LL[i+1]];
});
function segDist(a,b){ const dy=b[0]-a[0], dx=(b[1]-a[1])*0.758; return Math.hypot(dx,dy); }  // 0.758=cos(40.7N)
const segCum = segGeom.map(g=>{ const c=[0]; for(let k=1;k<g.length;k++) c.push(c[k-1]+segDist(g[k-1],g[k])); return c; });
function ptAlong(i,f){ const g=segGeom[i], c=segCum[i], tot=c[c.length-1];
  if(!tot || g.length<2) return g[0];
  const target=Math.max(0,Math.min(1,f))*tot; let lo=1,hi=c.length-1;
  while(lo<hi){ const m=(lo+hi)>>1; if(c[m]<target) lo=m+1; else hi=m; }
  const t0=c[lo-1], t1=c[lo], u=(t1>t0)?(target-t0)/(t1-t0):0;
  return [g[lo-1][0]+(g[lo][0]-g[lo-1][0])*u, g[lo-1][1]+(g[lo][1]-g[lo-1][1])*u]; }
function partialPts(i,f){ const g=segGeom[i], c=segCum[i], tot=c[c.length-1], p=ptAlong(i,f);
  if(!tot || g.length<2) return [g[0],p];
  const target=Math.max(0,Math.min(1,f))*tot, out=[];
  for(let k=0;k<g.length;k++){ if(c[k]<=target) out.push(g[k]); else break; }
  out.push(p); return out; }

// ridden track turns fully solid as the route passes (the faint network is the
// "not yet there" backdrop). train -> solid line color; run -> solid red dashed.
const doneSeg = [];
for (let i=0;i<WP.length-1;i++){
  const m = WP[i+1].mode;
  if (m==='train')
    doneSeg.push(L.polyline(segGeom[i], {renderer:routeR, color:segColor(i), weight:4, opacity:0,
      interactive:false, lineCap:'round', lineJoin:'round'}).addTo(map));
  else if (m==='run')
    doneSeg.push(L.polyline(segGeom[i], {renderer:routeR, color:'#e02b2b', weight:3, opacity:0,
      dashArray:'3 6', interactive:false, lineCap:'round'}).addTo(map));
  else doneSeg.push(null);                 // transfer / wait: no line
}
const partial = L.polyline([], {renderer:routeR, weight:4.5, opacity:1, lineCap:'round', lineJoin:'round'}).addTo(map);
const haloIcon = L.divIcon({className:'', html:'<div class="marker-halo" style="width:22px;height:22px"></div>', iconSize:[22,22]});
const halo = L.marker(LL[0], {icon:haloIcon, interactive:false, zIndexOffset:400}).addTo(map);
const marker = L.circleMarker(LL[0], {renderer:topR, radius:8, color:'#fff', weight:2.5, fillColor:'#111', fillOpacity:1}).addTo(map);

// --- visited-dot reconciliation (sorted by visit time; pointer in/out) ---
const visOrder = DATA.dots.map((d,idx)=>({idx,vis:d.vis})).filter(o=>o.vis!=null).sort((a,b)=>a.vis-b.vis);
let dp = 0;
function setVisited(e){
  while (dp<visOrder.length && visOrder[dp].vis<=e){ const o=visOrder[dp++]; dotLayers[o.idx].setStyle({...DOT_VIS, fillColor:DATA.dots[o.idx].vc||'#222'}); }
  while (dp>0 && visOrder[dp-1].vis>e){ const o=visOrder[--dp]; dotLayers[o.idx].setStyle(DOT_UNVIS); }
}

// --- readouts / legend / title ---
document.getElementById('total').textContent = DATA.total_hms;
document.getElementById('nstations').textContent = DATA.n_official;
document.getElementById('title').innerHTML =
  `<b>NYC Subway Challenge</b><br>${DATA.file} &middot; ${DATA.total_hms} &middot; ${DATA.n_official} stations`;
const legend = document.getElementById('legend');
legend.innerHTML = '<div style="font-weight:700;margin-bottom:4px">Lines ridden</div>';
for (const l of DATA.legend){ const s=document.createElement('span'); s.className='chip'; s.textContent=l.route;
  s.style.background=l.color; s.style.color=l.text; legend.append(s); }
const runLg=document.createElement('div'); runLg.style.marginTop='5px';
runLg.innerHTML='<span style="display:inline-block;width:16px;border-top:3px dashed #e02b2b;vertical-align:middle"></span> running';
legend.append(runLg);

function fmt(s){ s=Math.max(0,Math.round(s)); const h=(s/3600|0),m=(s%3600/60|0),x=s%60;
  return `${String(h).padStart(2,'0')}:${String(m).padStart(2,'0')}:${String(x).padStart(2,'0')}`; }
function clockOf(e){ const t=((DATA.start_t+e)%DATA.week+DATA.week)%DATA.week, day=DATA.days[(t/86400|0)%7], tod=t%86400;
  return [`${String(tod/3600|0).padStart(2,'0')}:${String(tod%3600/60|0).padStart(2,'0')}:${String(tod%60|0).padStart(2,'0')}`, day]; }
function modeLabel(b){ if(!b) return '—';
  if(b.mode==='train') return `● Riding the ${b.route}`;
  if(b.mode==='run') return '🏃 Running';
  if(b.mode==='transfer') return '⇄ Transfer';
  if(b.mode==='wait') return '⏸ Waiting'; return b.mode; }

// --- progress ---
const TOTAL = DATA.total_s;
let curIdx = 0;
function findSeg(e){ let lo=0,hi=WP.length-1; while(lo<hi){const m=(lo+hi+1)>>1; if(WP[m].e<=e) lo=m; else hi=m-1;} return Math.min(lo,WP.length-2); }
function setProgress(e){
  e=Math.max(0,Math.min(TOTAL,e));
  const i=findSeg(e), A=WP[i], B=WP[i+1];
  const span=Math.max(1,B.e-A.e), f=Math.max(0,Math.min(1,(e-A.e)/span));
  const pos=ptAlong(i,f);
  if (i>curIdx) for(let k=curIdx;k<i;k++){ if(doneSeg[k]) doneSeg[k].setStyle({opacity:1}); }
  else if (i<curIdx) for(let k=i;k<curIdx;k++){ if(doneSeg[k]) doneSeg[k].setStyle({opacity:0}); }
  curIdx=i;
  const isDwell = B.mode==='wait' || B.mode==='transfer';
  const mColor = B.mode==='run' ? '#e02b2b' : isDwell ? '#f5a623' : (B.color||'#111');
  if (B.mode==='train'){
    partial.setLatLngs(partialPts(i,f)); partial.setStyle({color:segColor(i), opacity:1, dashArray:null, weight:4.5});
  } else if (B.mode==='run'){
    partial.setLatLngs(partialPts(i,f)); partial.setStyle({color:'#e02b2b', opacity:1, dashArray:'3 6', weight:3});
  } else partial.setLatLngs([]);
  marker.setLatLng(pos).setStyle({fillColor:mColor});
  halo.setLatLng(pos);
  const hd=halo.getElement(); if(hd&&hd.firstElementChild) hd.firstElementChild.style.setProperty('--pc', rgba(mColor,0.6));
  setVisited(e);
  const [clk,day]=clockOf(e);
  document.getElementById('clock').textContent=clk;
  document.getElementById('day').textContent=day;
  document.getElementById('elapsed').textContent=fmt(e);
  document.getElementById('visited').textContent=dp;
  const st=document.getElementById('status'); st.textContent=modeLabel(B) + (isDwell ? ` ${fmtDur(B.e-A.e)}` : '');
  st.style.color = mColor;
  slider.value=Math.round(e/TOTAL*1000);
  updateBreakdown(e);
  syncItin(e);
}

// --- controls ---
const slider=document.getElementById('slider'), playBtn=document.getElementById('play');
let e=0, playing=false, mult=800, raf=0, last=0;
slider.addEventListener('input',()=>{ e=slider.value/1000*TOTAL; setProgress(e); });
function tick(ts){ if(!playing) return;
  if(last){ e+=(ts-last)/1000*mult; if(e>=TOTAL){ e=TOTAL; playing=false; playBtn.textContent='▶ Play'; playBtn.classList.remove('on'); } }
  last=ts; setProgress(e); if(playing) raf=requestAnimationFrame(tick); }
playBtn.addEventListener('click',()=>{ playing=!playing; playBtn.textContent=playing?'⏸ Pause':'▶ Play'; playBtn.classList.toggle('on',playing);
  if(playing){ if(e>=TOTAL) e=0; last=0; raf=requestAnimationFrame(tick); } });
for(const b of document.querySelectorAll('.speed')) b.addEventListener('click',()=>{ mult=+b.dataset.mult;
  document.querySelectorAll('.speed').forEach(x=>x.classList.remove('on')); b.classList.add('on'); });

// --- itinerary side panel (synced trip log) ---
function txtOn(hex){ const c=hex.replace('#',''); const r=parseInt(c.substr(0,2),16),g=parseInt(c.substr(2,2),16),b=parseInt(c.substr(4,2),16);
  return (0.299*r+0.587*g+0.114*b)>150 ? '#000' : '#fff'; }
function fmtDur(s){ s=Math.round(s); if(s>=3600){const h=s/3600|0,m=(s%3600/60|0); return `${h}h${String(m).padStart(2,'0')}`;} if(s>=60) return `${(s/60|0)}m`; return `${s}s`; }
const itList=document.getElementById('itlist'), itE0=DATA.itinerary.map(o=>o.e0);
document.getElementById('itcount').textContent = `· ${DATA.itinerary.length} steps`;
const itEls = DATA.itinerary.map(o=>{
  const c0=clockOf(o.e0)[0], dur=fmtDur(o.e1-o.e0);
  let badge,title,sub;
  if(o.type==='ride'){ badge=o.line; title=`Ride the ${o.line}`; sub=`${c0} · ${dur} · ${o.stops} stop${o.stops>1?'s':''} → ${o.to}`; }
  else if(o.type==='run'){ badge='🏃'; title='Run'; sub=`${c0} · ${dur} · ${o.frm} → ${o.to}`; }
  else { badge='⇄'; title=o.label; sub=`${c0} · ${dur} · ${o.to}`; }
  const div=document.createElement('div'); div.className='it';
  div.innerHTML=`<div class="ic" style="background:${o.color};color:${txtOn(o.color)}">${badge}</div>`+
                `<div><div class="t">${title}</div><div class="s">${sub}</div></div>`;
  div.addEventListener('click',()=>{ e=o.e0; setProgress(e); });
  itList.append(div); return div;
});
let curIt=-1;
function syncItin(ev){
  let lo=0,hi=itE0.length; while(lo<hi){const m=(lo+hi)>>1; if(itE0[m]<=ev) lo=m+1; else hi=m;}
  const i=Math.max(0,lo-1);
  if(i!==curIt){ if(curIt>=0) itEls[curIt].classList.remove('cur');
    itEls[i].classList.add('cur'); itEls[i].scrollIntoView({block:'nearest',behavior:'smooth'}); curIt=i; }
}

setTimeout(()=>map.invalidateSize(), 100);
setProgress(0);
</script>
</body>
</html>
"""


def render_html(payload: dict) -> str:
    return _LEAFLET_HTML.replace("/*__DATA__*/", json.dumps(payload, separators=(",", ":")))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Render a solution as an interactive HTML replay.")
    p.add_argument("solution", nargs="?", default="solutions/best.json",
                   help="Solution JSON (default: solutions/best.json)")
    p.add_argument("-o", "--out", default=None,
                   help="Output .html (default: alongside the solution, _replay_map.html)")
    p.add_argument("--radius", type=int, default=None, help="Run-layer radius (m).")
    args = p.parse_args(argv)

    src = Path(args.solution)
    out = Path(args.out) if args.out else src.with_name(src.stem + "_replay_map.html")
    payload = build_payload(src, radius_m=args.radius)
    out.write_text(render_html(payload))
    print(f"wrote {out}  ({out.stat().st_size//1024} KB)")
    print(f"  route: {payload['total_hms']}  stations {payload['n_visited']}/{payload['n_official']}"
          f"  waypoints {len(payload['waypoints'])}  network lines {len(payload['shapes'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Small NEOS XML-RPC client for async optimization experiments.

The goal is not to hide NEOS behind a large framework. It is to make compact
solver experiments reproducible: submit an XML job, record its job number and
password in ignored artifacts, poll it later, and fetch results without printing
secrets from ``.env``.
"""
from __future__ import annotations

import argparse
import json
import os
import ssl
import sys
import time
import xmlrpc.client
from datetime import datetime, timezone
from pathlib import Path

import certifi

NEOS_URL = "https://neos-server.org:3333"
DEFAULT_ARTIFACT_DIR = Path("reports/optimization_runs")
JOBS_FILE = "neos_jobs.jsonl"


def load_dotenv(path: Path = Path(".env")) -> dict[str, str]:
    """Load simple KEY=VALUE entries from .env without echoing values."""
    loaded = {}
    if not path.exists():
        return loaded
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip("'\"")
        if key and key not in os.environ:
            os.environ[key] = value
            loaded[key] = value
    return loaded


def artifact_dir() -> Path:
    out = Path(os.environ.get("OPT_ARTIFACT_DIR", DEFAULT_ARTIFACT_DIR))
    out.mkdir(parents=True, exist_ok=True)
    return out


def neos_server(url: str = NEOS_URL):
    context = ssl.create_default_context(cafile=certifi.where())
    transport = xmlrpc.client.SafeTransport(context=context)
    return xmlrpc.client.ServerProxy(url, transport=transport, allow_none=True)


def _decode_blob(blob) -> str:
    if hasattr(blob, "data"):
        data = blob.data
    else:
        data = blob
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    return str(data)


def _jobs_path() -> Path:
    return artifact_dir() / JOBS_FILE


def record_job(job_number: int, password: str, label: str, source: str) -> Path:
    entry = {
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "job_number": int(job_number),
        "password": str(password),
        "label": label,
        "source": source,
    }
    path = _jobs_path()
    with path.open("a") as f:
        f.write(json.dumps(entry, sort_keys=True) + "\n")
    return path


def latest_job() -> dict | None:
    path = _jobs_path()
    if not path.exists():
        return None
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    return rows[-1] if rows else None


def _job_args(args):
    if args.latest:
        row = latest_job()
        if not row:
            raise SystemExit(f"no recorded NEOS jobs in {_jobs_path()}")
        return int(row["job_number"]), str(row["password"])
    if args.job is None or args.password is None:
        raise SystemExit("provide --job and --password, or use --latest")
    return int(args.job), str(args.password)


def ampl_xml(model: str, data: str, commands: str, email: str,
             solver: str = "Cbc", category: str = "milp") -> str:
    return f"""<document>
<category>{category}</category>
<solver>{solver}</solver>
<inputMethod>AMPL</inputMethod>
<email>{email}</email>

<model><![CDATA[{model}]]></model>

<data><![CDATA[{data}]]></data>

<commands><![CDATA[{commands}]]></commands>

<comments><![CDATA[Submitted by subway_challenge.neos_client]]></comments>

</document>
"""


def tiny_ampl_job(email: str) -> str:
    model = """
var x binary;
maximize obj: x;
subject to limit: x <= 1;
"""
    commands = """
option solver cbc;
solve;
display x, obj;
"""
    return ampl_xml(model=model, data="", commands=commands, email=email)


def cmd_ping(args) -> int:
    load_dotenv()
    server = neos_server(args.url)
    print(server.ping())
    if args.solvers:
        solvers = list(server.listAllSolvers())
        print(f"solvers={len(solvers)}")
        for item in solvers[:args.solvers]:
            print(item)
    return 0


def cmd_template(args) -> int:
    load_dotenv()
    server = neos_server(args.url)
    text = server.getSolverTemplate(args.category, args.solver, args.input_method)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
        print(f"wrote {out}")
    else:
        print(text)
    return 0


def submit_xml(server, xml_text: str, label: str, source: str, wait: bool) -> tuple[int, str]:
    job_number, password = server.submitJob(xml_text)
    try:
        job_number = int(job_number)
    except (TypeError, ValueError):
        job_number = 0
    if job_number <= 0:
        raise SystemExit(f"NEOS submission failed: {password}")
    path = record_job(job_number, str(password), label, source)
    print(f"submitted job={job_number} label={label!r}; recorded in {path}")
    if wait:
        while True:
            status = server.getJobStatus(job_number, password)
            print(f"status={status}", flush=True)
            if status == "Done":
                break
            time.sleep(10)
        result = _decode_blob(server.getFinalResults(job_number, password))
        out = artifact_dir() / f"neos_{job_number}_{label}.txt"
        out.write_text(result)
        print(f"wrote {out}")
    return job_number, str(password)


def cmd_submit_xml(args) -> int:
    load_dotenv()
    xml_path = Path(args.xml)
    server = neos_server(args.url)
    submit_xml(
        server,
        xml_path.read_text(),
        label=args.label or xml_path.stem,
        source=str(xml_path),
        wait=args.wait,
    )
    return 0


def cmd_submit_tiny(args) -> int:
    load_dotenv()
    email = os.environ.get("NEOS_EMAIL")
    if not email:
        raise SystemExit("NEOS_EMAIL must be set in .env or the environment")
    server = neos_server(args.url)
    submit_xml(
        server,
        tiny_ampl_job(email),
        label=args.label,
        source="tiny-ampl",
        wait=args.wait,
    )
    return 0


def cmd_status(args) -> int:
    load_dotenv()
    job, password = _job_args(args)
    server = neos_server(args.url)
    status = server.getJobStatus(job, password)
    print(f"job={job} status={status}")
    if args.info:
        print(server.getJobInfo(job, password))
    return 0


def cmd_fetch(args) -> int:
    load_dotenv()
    job, password = _job_args(args)
    server = neos_server(args.url)
    if args.blocking:
        text = _decode_blob(server.getFinalResults(job, password))
    else:
        text = _decode_blob(server.getFinalResultsNonBlocking(job, password))
        if not text:
            print(f"job={job} results not ready")
            return 2
    label = args.label or f"job_{job}"
    out = Path(args.out) if args.out else artifact_dir() / f"neos_{job}_{label}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(f"wrote {out}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="NEOS XML-RPC helper for optimization runs.")
    p.add_argument("--url", default=NEOS_URL)
    sub = p.add_subparsers(dest="cmd", required=True)

    ping = sub.add_parser("ping", help="Ping NEOS and optionally list solvers.")
    ping.add_argument("--solvers", type=int, default=0,
                      help="Print the first N solver interfaces after pinging.")
    ping.set_defaults(func=cmd_ping)

    tmpl = sub.add_parser("template", help="Fetch a NEOS solver XML template.")
    tmpl.add_argument("category")
    tmpl.add_argument("solver")
    tmpl.add_argument("input_method")
    tmpl.add_argument("--out")
    tmpl.set_defaults(func=cmd_template)

    sx = sub.add_parser("submit-xml", help="Submit a prepared NEOS XML job.")
    sx.add_argument("xml")
    sx.add_argument("--label")
    sx.add_argument("--wait", action="store_true")
    sx.set_defaults(func=cmd_submit_xml)

    tiny = sub.add_parser("submit-tiny-ampl", help="Submit a tiny Cbc/AMPL smoke job.")
    tiny.add_argument("--label", default="tiny_ampl")
    tiny.add_argument("--wait", action="store_true")
    tiny.set_defaults(func=cmd_submit_tiny)

    status = sub.add_parser("status", help="Check a NEOS job status.")
    status.add_argument("--job", type=int)
    status.add_argument("--password")
    status.add_argument("--latest", action="store_true")
    status.add_argument("--info", action="store_true")
    status.set_defaults(func=cmd_status)

    fetch = sub.add_parser("fetch", help="Fetch final NEOS results.")
    fetch.add_argument("--job", type=int)
    fetch.add_argument("--password")
    fetch.add_argument("--latest", action="store_true")
    fetch.add_argument("--blocking", action="store_true")
    fetch.add_argument("--label")
    fetch.add_argument("--out")
    fetch.set_defaults(func=cmd_fetch)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

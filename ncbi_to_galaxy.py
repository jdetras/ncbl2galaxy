#!/usr/bin/env python3
import argparse
import json
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import requests
from defusedxml import ElementTree as ET
from requests import RequestException
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
RUN_REGEX = re.compile(r"\b([SED]RR\d+)\b")


@dataclass
class RunRecord:
    run_accession: str
    sample_accession: str
    library_layout: str
    fastq_urls: list[str]


class StateStore:
    def __init__(self, path: str):
        self.path = Path(path)
        self.processed_runs: set[str] = set()

    def load(self) -> None:
        if not self.path.exists():
            return
        data = json.loads(self.path.read_text(encoding="utf-8"))
        self.processed_runs = set(data.get("processed_runs", []))

    def save(self) -> None:
        payload = {"processed_runs": sorted(self.processed_runs)}
        self.path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def build_retry_session(total_retries: int = 5, backoff_factor: float = 1.0) -> requests.Session:
    retry = Retry(
        total=total_retries,
        read=total_retries,
        connect=total_retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def chunked(values: list[str], size: int) -> list[list[str]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def ncbi_get(
    session: requests.Session,
    endpoint: str,
    params: dict[str, str],
    timeout: int = 60,
) -> requests.Response:
    resp = session.get(f"{EUTILS_BASE}/{endpoint}", params=params, timeout=timeout)
    resp.raise_for_status()
    return resp


def esearch_pubmed(
    session: requests.Session,
    query: str,
    retmax: int,
    email: str | None,
    ncbi_api_key: str | None,
) -> list[str]:
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": str(retmax),
        "retmode": "json",
    }
    if email:
        params["email"] = email
    if ncbi_api_key:
        params["api_key"] = ncbi_api_key

    data = ncbi_get(session, "esearch.fcgi", params).json()
    return data.get("esearchresult", {}).get("idlist", [])


def elink_pubmed_to_sra(
    session: requests.Session,
    pmids: list[str],
    email: str | None,
    ncbi_api_key: str | None,
) -> set[str]:
    sra_ids: set[str] = set()
    for batch in chunked(pmids, 200):
        params = {
            "dbfrom": "pubmed",
            "db": "sra",
            "id": ",".join(batch),
            "retmode": "xml",
            "linkname": "pubmed_sra",
        }
        if email:
            params["email"] = email
        if ncbi_api_key:
            params["api_key"] = ncbi_api_key

        root = ET.fromstring(ncbi_get(session, "elink.fcgi", params).text)
        for node in root.findall(".//LinkSetDb/Link/Id"):
            if node.text:
                sra_ids.add(node.text.strip())
        time.sleep(0.34)
    return sra_ids


def esummary_sra_runs(
    session: requests.Session,
    sra_ids: list[str],
    email: str | None,
    ncbi_api_key: str | None,
) -> set[str]:
    runs: set[str] = set()
    for batch in chunked(sra_ids, 200):
        params = {
            "db": "sra",
            "id": ",".join(batch),
            "retmode": "xml",
        }
        if email:
            params["email"] = email
        if ncbi_api_key:
            params["api_key"] = ncbi_api_key

        root = ET.fromstring(ncbi_get(session, "esummary.fcgi", params).text)
        for item in root.findall(".//DocSum/Item[@Name='Runs']"):
            text = item.text or ""
            for match in RUN_REGEX.findall(text):
                runs.add(match)
        time.sleep(0.34)
    return runs


def ena_run_record(
    session: requests.Session, run_accession: str, timeout: int = 60
) -> RunRecord | None:
    params = {
        "accession": run_accession,
        "result": "read_run",
        "fields": "fastq_ftp,sample_accession,library_layout",
        "format": "json",
    }
    resp = session.get(ENA_FILEREPORT, params=params, timeout=timeout)
    resp.raise_for_status()

    rows = resp.json()
    if not rows:
        return None

    row = rows[0]
    fastq_ftp = row.get("fastq_ftp", "")
    if not fastq_ftp:
        return None

    urls: list[str] = []
    for part in fastq_ftp.split(";"):
        part = part.strip()
        if not part:
            continue
        if part.startswith("http://") or part.startswith("https://"):
            urls.append(part)
        elif part.startswith("ftp."):
            urls.append(f"https://{part}")
        else:
            urls.append(f"https://{part}")

    if not urls:
        return None

    sample_accession = row.get("sample_accession") or run_accession
    library_layout = (row.get("library_layout") or "").upper() or (
        "PAIRED" if len(urls) >= 2 else "SINGLE"
    )

    return RunRecord(
        run_accession=run_accession,
        sample_accession=sample_accession,
        library_layout=library_layout,
        fastq_urls=urls,
    )


class GalaxyClient:
    def __init__(self, base_url: str, api_key: str, session: requests.Session | None = None):
        self.base_url = base_url.rstrip("/")
        self.session = session or build_retry_session()
        self.session.headers.update({"x-api-key": api_key, "Content-Type": "application/json"})

    def _get(self, path: str, params: dict[str, str] | None = None):
        resp = self.session.get(f"{self.base_url}{path}", params=params, timeout=120)
        resp.raise_for_status()
        return resp.json()

    def _post(self, path: str, payload: dict):
        resp = self.session.post(f"{self.base_url}{path}", data=json.dumps(payload), timeout=120)
        resp.raise_for_status()
        return resp.json()

    def list_workflows(self) -> list[dict]:
        return self._get("/api/workflows")

    def find_workflow_id_by_name(self, workflow_name: str) -> str:
        workflows = self.list_workflows()
        matches = [wf for wf in workflows if wf.get("name") == workflow_name]
        if not matches:
            raise RuntimeError(f"Workflow named '{workflow_name}' was not found")
        if len(matches) > 1:
            raise RuntimeError(
                f"Multiple workflows named '{workflow_name}' found; use an explicit workflow ID"
            )
        return matches[0]["id"]

    def create_history(self, name: str) -> str:
        data = self._post("/api/histories", {"name": name})
        return data["id"]

    def get_workflow_input_id(self, workflow_id: str, input_label: str | None) -> str:
        wf = self._get(f"/api/workflows/{workflow_id}")
        inputs = wf.get("inputs", {})
        if not inputs:
            raise RuntimeError("Workflow has no declared inputs")

        if input_label:
            for input_id, meta in inputs.items():
                if meta.get("label") == input_label or meta.get("name") == input_label:
                    return input_id
            raise RuntimeError(f"Workflow input label '{input_label}' not found")

        return sorted(inputs.keys(), key=lambda x: int(x) if x.isdigit() else x)[0]

    def fetch_url_to_history(self, history_id: str, url: str, name: str) -> str:
        payload = {
            "history_id": history_id,
            "targets": [
                {
                    "destination": {"type": "hdas"},
                    "elements": [{"src": "url", "url": url, "name": name}],
                }
            ],
        }
        data = self._post("/api/tools/fetch", payload)
        outputs = data.get("outputs", [])
        if not outputs:
            raise RuntimeError(f"No outputs returned when fetching {url}")
        return outputs[0]["id"]

    def create_pair_collection(
        self, history_id: str, forward_id: str, reverse_id: str, name: str
    ) -> str:
        payload = {
            "history_id": history_id,
            "collection_type": "paired",
            "name": name,
            "element_identifiers": [
                {"name": "forward", "src": "hda", "id": forward_id},
                {"name": "reverse", "src": "hda", "id": reverse_id},
            ],
        }
        data = self._post("/api/dataset_collections", payload)
        return data["id"]

    def invoke_workflow(
        self, workflow_id: str, history_id: str, inputs: dict[str, dict[str, str]]
    ) -> str:
        payload = {
            "history_id": history_id,
            "inputs": inputs,
        }
        data = self._post(f"/api/workflows/{workflow_id}/invocations", payload)
        return data.get("id", "")


def group_runs_by_sample(run_records: list[RunRecord]) -> dict[str, list[RunRecord]]:
    grouped: dict[str, list[RunRecord]] = {}
    for record in run_records:
        grouped.setdefault(record.sample_accession, []).append(record)
    return grouped


def resolve_workflow_id(
    client: GalaxyClient, workflow_id: str | None, workflow_name: str | None
) -> str | None:
    if workflow_id:
        return workflow_id
    if workflow_name:
        return client.find_workflow_id_by_name(workflow_name)
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Search NCBI PubMed for rice-related publications, discover linked SRA runs, "
            "and optionally import reads into Galaxy and invoke variant-calling workflows."
        )
    )
    parser.add_argument("--query", default='(rice[Title/Abstract] OR "Oryza sativa"[MeSH Terms])')
    parser.add_argument("--retmax", type=int, default=50)
    parser.add_argument("--max-runs", type=int, default=20)
    parser.add_argument("--email", help="Email passed to NCBI E-utilities")
    parser.add_argument("--ncbi-api-key", help="NCBI API key for higher rate limits")

    parser.add_argument("--galaxy-url", help="Galaxy base URL, e.g. https://usegalaxy.org")
    parser.add_argument("--galaxy-api-key", help="Galaxy user API key")

    parser.add_argument("--single-workflow-id", help="Galaxy single-end workflow ID")
    parser.add_argument(
        "--single-workflow-name",
        default="Rice Variant Calling (BWA-MEM2 + FreeBayes)",
        help="Galaxy single-end workflow name if ID not supplied",
    )
    parser.add_argument("--paired-workflow-id", help="Galaxy paired-end workflow ID")
    parser.add_argument(
        "--paired-workflow-name",
        default="Rice Variant Calling Paired (BWA-MEM2 + FreeBayes)",
        help="Galaxy paired-end workflow name if ID not supplied",
    )

    parser.add_argument(
        "--single-input-label", default="Reads FASTQ", help="Single-end workflow read input"
    )
    parser.add_argument(
        "--paired-input-label", default="Reads Pair", help="Paired-end workflow read input"
    )
    parser.add_argument(
        "--reference-input-label",
        default="Reference FASTA",
        help="Workflow label/name for reference input",
    )
    parser.add_argument("--reference-url", help="Optional reference FASTA URL to fetch into Galaxy")
    parser.add_argument(
        "--reference-dataset-id", help="Optional existing Galaxy dataset ID for reference FASTA"
    )

    parser.add_argument("--history-id", help="Existing Galaxy history ID")
    parser.add_argument("--history-name", default="Rice_variant_calling_inputs")
    parser.add_argument(
        "--history-per-sample",
        action="store_true",
        help="Create one Galaxy history per sample accession",
    )

    parser.add_argument(
        "--state-file", default=".ncbi_to_galaxy_state.json", help="State file for resume mode"
    )
    parser.add_argument(
        "--reset-state", action="store_true", help="Ignore previous state and start fresh"
    )
    parser.add_argument("--dry-run", action="store_true", help="Do not upload or invoke workflow")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    state = StateStore(args.state_file)
    if not args.reset_state:
        state.load()

    ncbi_session = build_retry_session()
    print(f"Searching PubMed with query: {args.query}", flush=True)
    try:
        pmids = esearch_pubmed(ncbi_session, args.query, args.retmax, args.email, args.ncbi_api_key)
    except RequestException as exc:
        print(f"ERROR failed to query NCBI PubMed: {exc}", file=sys.stderr)
        return 1

    print(f"Found {len(pmids)} PubMed records")
    if not pmids:
        return 0

    try:
        sra_ids = sorted(elink_pubmed_to_sra(ncbi_session, pmids, args.email, args.ncbi_api_key))
    except RequestException as exc:
        print(f"ERROR failed to map PubMed records to SRA: {exc}", file=sys.stderr)
        return 1

    print(f"Found {len(sra_ids)} linked SRA record IDs")
    if not sra_ids:
        return 0

    try:
        run_accessions = sorted(
            esummary_sra_runs(ncbi_session, sra_ids, args.email, args.ncbi_api_key)
        )
    except RequestException as exc:
        print(f"ERROR failed to resolve SRA run accessions: {exc}", file=sys.stderr)
        return 1

    if args.max_runs > 0:
        run_accessions = run_accessions[: args.max_runs]

    print(f"Resolved {len(run_accessions)} run accessions (capped by --max-runs)")
    if not run_accessions:
        return 0

    run_records: list[RunRecord] = []
    for run in run_accessions:
        if run in state.processed_runs:
            continue
        try:
            record = ena_run_record(ncbi_session, run)
            if record:
                run_records.append(record)
        except RequestException as exc:
            print(f"WARN failed to resolve FASTQ URLs for {run}: {exc}", file=sys.stderr)

    print(f"Runs with FASTQ URLs and not already processed: {len(run_records)}")
    if not run_records:
        return 0

    grouped = group_runs_by_sample(run_records)
    print(f"Grouped into {len(grouped)} biological samples")

    for sample, records in list(grouped.items())[:10]:
        layouts = {r.library_layout for r in records}
        print(f"  {sample}: {len(records)} run(s), layouts={','.join(sorted(layouts))}")

    if args.dry_run:
        return 0

    galaxy_ready = all([args.galaxy_url, args.galaxy_api_key])
    if not galaxy_ready:
        print("Galaxy parameters missing; provide --galaxy-url and --galaxy-api-key.")
        return 0

    galaxy = GalaxyClient(args.galaxy_url, args.galaxy_api_key, session=build_retry_session())

    single_workflow_id = resolve_workflow_id(
        galaxy, args.single_workflow_id, args.single_workflow_name
    )
    paired_workflow_id = resolve_workflow_id(
        galaxy, args.paired_workflow_id, args.paired_workflow_name
    )

    if not single_workflow_id:
        print("ERROR no single-end workflow configured", file=sys.stderr)
        return 1

    single_input_id = galaxy.get_workflow_input_id(single_workflow_id, args.single_input_label)
    reference_input_single_id = galaxy.get_workflow_input_id(
        single_workflow_id, args.reference_input_label
    )

    paired_input_id = None
    reference_input_paired_id = None
    if paired_workflow_id:
        paired_input_id = galaxy.get_workflow_input_id(paired_workflow_id, args.paired_input_label)
        reference_input_paired_id = galaxy.get_workflow_input_id(
            paired_workflow_id, args.reference_input_label
        )

    history_cache: dict[str, str] = {}

    def history_for_sample(sample_accession: str) -> str:
        if args.history_id:
            return args.history_id
        if not args.history_per_sample:
            key = "shared"
            if key not in history_cache:
                history_cache[key] = galaxy.create_history(args.history_name)
            return history_cache[key]
        key = sample_accession
        if key not in history_cache:
            history_cache[key] = galaxy.create_history(f"{args.history_name}_{sample_accession}")
        return history_cache[key]

    reference_cache: dict[str, str] = {}

    def reference_for_history(history_id: str) -> str | None:
        if args.reference_dataset_id:
            return args.reference_dataset_id
        if not args.reference_url:
            return None
        if history_id not in reference_cache:
            reference_cache[history_id] = galaxy.fetch_url_to_history(
                history_id, args.reference_url, "reference.fasta"
            )
        return reference_cache[history_id]

    uploaded = 0
    invoked = 0

    for sample_accession, records in grouped.items():
        history_id = history_for_sample(sample_accession)
        ref_dataset_id = reference_for_history(history_id)

        for record in records:
            try:
                if len(record.fastq_urls) >= 2 or record.library_layout == "PAIRED":
                    if not (paired_workflow_id and paired_input_id and reference_input_paired_id):
                        print(
                            f"WARN skipping paired run {record.run_accession}: paired workflow is not configured",
                            file=sys.stderr,
                        )
                        continue

                    fwd = galaxy.fetch_url_to_history(
                        history_id,
                        record.fastq_urls[0],
                        f"{record.run_accession}_R1.fastq.gz",
                    )
                    rev = galaxy.fetch_url_to_history(
                        history_id,
                        record.fastq_urls[1],
                        f"{record.run_accession}_R2.fastq.gz",
                    )
                    pair_id = galaxy.create_pair_collection(
                        history_id, fwd, rev, f"{record.run_accession}_pair"
                    )
                    uploaded += 2

                    workflow_inputs = {paired_input_id: {"src": "hdca", "id": pair_id}}
                    if ref_dataset_id:
                        workflow_inputs[reference_input_paired_id] = {
                            "src": "hda",
                            "id": ref_dataset_id,
                        }

                    invocation_id = galaxy.invoke_workflow(
                        paired_workflow_id, history_id, workflow_inputs
                    )
                    invoked += 1
                    print(
                        f"Sample {sample_accession}: paired run {record.run_accession} uploaded, invocation={invocation_id}"
                    )
                else:
                    read_id = galaxy.fetch_url_to_history(
                        history_id,
                        record.fastq_urls[0],
                        f"{record.run_accession}.fastq.gz",
                    )
                    uploaded += 1
                    workflow_inputs = {single_input_id: {"src": "hda", "id": read_id}}
                    if ref_dataset_id:
                        workflow_inputs[reference_input_single_id] = {
                            "src": "hda",
                            "id": ref_dataset_id,
                        }

                    invocation_id = galaxy.invoke_workflow(
                        single_workflow_id, history_id, workflow_inputs
                    )
                    invoked += 1
                    print(
                        f"Sample {sample_accession}: single run {record.run_accession} uploaded, invocation={invocation_id}"
                    )

                state.processed_runs.add(record.run_accession)
                state.save()
            except Exception as exc:
                print(f"ERROR processing run {record.run_accession}: {exc}", file=sys.stderr)

    print(f"Completed. Uploaded {uploaded} datasets and invoked {invoked} workflows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

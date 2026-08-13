#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_ROOT="${1:-${REPO_ROOT}/results/orchestration-smoke/$(date +%Y%m%d-%H%M%S)}"

if [[ -e "${RUN_ROOT}" ]]; then
    echo "Refusing to overwrite existing smoke directory: ${RUN_ROOT}" >&2
    exit 64
fi

if [[ -n "${NEXTFLOW_BIN:-}" ]]; then
    NEXTFLOW="${NEXTFLOW_BIN}"
elif command -v nextflow >/dev/null 2>&1; then
    NEXTFLOW="$(command -v nextflow)"
elif [[ -x "${REPO_ROOT}/.tools/nextflow" ]]; then
    NEXTFLOW="${REPO_ROOT}/.tools/nextflow"
else
    echo "Nextflow is required for the orchestration smoke test." >&2
    exit 69
fi

if [[ -d "${REPO_ROOT}/.tools/nextflow-env" ]]; then
    export JAVA_HOME="${JAVA_HOME:-${REPO_ROOT}/.tools/nextflow-env}"
fi
export NXF_HOME="${NXF_HOME:-${REPO_ROOT}/.tools/nextflow-home}"
export NXF_VER="${NXF_VER:-24.10.4}"
export NXF_OFFLINE="${NXF_OFFLINE:-true}"

mkdir -p "${RUN_ROOT}"

run_batch() {
    local batch_number="$1"
    local query_file="$2"
    local batch_root="${RUN_ROOT}/batch-${batch_number}"
    local delivery="${batch_root}/delivery"
    local work="${batch_root}/work"

    mkdir -p "${batch_root}/audit"

    local command=(
        "${NEXTFLOW}"
        run "${REPO_ROOT}/workflow/main.nf"
        -profile test
        --queries "${query_file}"
        --outdir "${delivery}"
        -work-dir "${work}"
    )

    "${command[@]}" \
        -with-report "${batch_root}/audit/first-report.html" \
        -with-trace "${batch_root}/audit/first-trace.tsv" \
        -with-timeline "${batch_root}/audit/first-timeline.html" \
        2>&1 | tee "${batch_root}/first-run.log"

    "${command[@]}" -resume \
        -with-report "${batch_root}/audit/resume-report.html" \
        -with-trace "${batch_root}/audit/resume-trace.tsv" \
        -with-timeline "${batch_root}/audit/resume-timeline.html" \
        2>&1 | tee "${batch_root}/resume-run.log"
}

run_batch 1 "${REPO_ROOT}/data/benchmark/smoke_v1_batch_a.smi"
run_batch 2 "${REPO_ROOT}/data/benchmark/smoke_v1_batch_b.smi"

python - "${RUN_ROOT}" <<'PY'
import hashlib
import json
import pathlib
import sys

run_root = pathlib.Path(sys.argv[1])
manifests = []
resume_logs = []

for batch_number in (1, 2):
    batch_root = run_root / f"batch-{batch_number}"
    manifest_path = batch_root / "delivery/final_report/report_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    report_path = manifest_path.parent / manifest["report"]["path"]
    report_digest = hashlib.sha256(report_path.read_bytes()).hexdigest()
    if report_digest != manifest["report"]["sha256"]:
        raise SystemExit(f"report hash mismatch: {report_path}")
    if not 1 <= int(manifest["query_count"]) <= 5:
        raise SystemExit(f"invalid service batch size: {manifest['query_count']}")
    if manifest["validation_scope"] != "orchestration_mock_only":
        raise SystemExit("mock smoke produced an unexpected validation scope")
    manifests.append(manifest)
    resume_logs.append((batch_root / "resume-run.log").read_text(encoding="utf-8"))

queries = [query for manifest in manifests for query in manifest["queries"]]
successes = sum(
    int(manifest["query_count"]) * float(manifest["technical_success_rate"])
    for manifest in manifests
)
resumed_without_recomputation = all(
    "cached:" in log.lower() and "submitted process" not in log.lower()
    for log in resume_logs
)

summary = {
    "schema_version": "1.0",
    "validation_scope": "orchestration_mock_only",
    "batch_count": len(manifests),
    "query_count": len(queries),
    "max_queries_per_batch": max(int(item["query_count"]) for item in manifests),
    "technical_success_rate": successes / len(queries),
    "all_metrics_traceable": all(item["all_metrics_traceable"] for item in manifests),
    "unsupported_targets_have_no_numeric_score": all(
        item["unsupported_targets_have_no_numeric_score"] for item in manifests
    ),
    "resumed_without_recomputation": resumed_without_recomputation,
    "queries_are_unique": len(queries) == len(set(queries)),
    "queries": queries,
}

if summary["query_count"] != 10 or not summary["queries_are_unique"]:
    raise SystemExit("the orchestration smoke must contain ten unique queries")
if summary["technical_success_rate"] < 0.95:
    raise SystemExit("technical success rate is below 95%")
if not summary["all_metrics_traceable"]:
    raise SystemExit("a mock report metric is not traceable")
if not summary["unsupported_targets_have_no_numeric_score"]:
    raise SystemExit("unsupported-target missingness was not preserved")
if not summary["resumed_without_recomputation"]:
    raise SystemExit("Nextflow resume recomputed at least one process")

(run_root / "orchestration_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY

(
    cd "${RUN_ROOT}"
    find . -type f ! -name artifacts.sha256 -print0 \
        | sort -z \
        | xargs -0 sha256sum > artifacts.sha256
)

echo "AIRTI orchestration-only smoke passed: ${RUN_ROOT}"

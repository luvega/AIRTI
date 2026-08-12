"""Render release-gated Markdown reports with artifact-level provenance."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, StrictUndefined, select_autoescape
from pydantic import BaseModel, ConfigDict, Field, model_validator

from airti_tf.manifest_io import content_sha256, write_bytes_atomic

TEMPLATE_ROOT = Path(__file__).resolve().parents[3] / "templates"
PROHIBITED_CLAIMS = ("确认靶点", "已证实直接结合", "真实结合概率")


class ProhibitedClaimError(ValueError):
    """Raised when a computation-only release overstates its evidence."""


class ArtifactIntegrityError(ValueError):
    """Raised when a traced artifact is missing or has a mismatched digest."""


class ArtifactLink(BaseModel):
    model_config = ConfigDict(extra="forbid")

    artifact_id: str
    task_id: str
    path: Path
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class CoverageContext(BaseModel):
    model_config = ConfigDict(extra="forbid")

    total: int = Field(ge=0)
    ready: int = Field(ge=0)
    unsupported: int = Field(ge=0)
    artifact_id: str


class TargetContext(BaseModel):
    model_config = ConfigDict(extra="forbid")

    rank: int = Field(gt=0)
    target_id: str
    gene_symbol: str
    priority: float = Field(ge=0, le=1)
    evidence_tier: str
    vina_score: float | None = None
    boltz_score: float | None = None
    md_score: float | None = None
    uncertainty_flags: list[str] = Field(default_factory=list)
    artifact_id: str


class ReportContext(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: str
    generated_at: str
    query_ligands: list[str] = Field(min_length=1, max_length=5)
    conclusion: str
    coverage: CoverageContext
    top_targets: list[TargetContext]
    artifacts: list[ArtifactLink]
    limitations: list[str]
    wet_lab_recommendations: list[str]

    @model_validator(mode="after")
    def validate_artifact_references(self) -> "ReportContext":
        declared = {artifact.artifact_id for artifact in self.artifacts}
        if len(declared) != len(self.artifacts):
            raise ValueError("artifact_id values must be unique")
        referenced = {self.coverage.artifact_id} | {
            target.artifact_id for target in self.top_targets
        }
        missing = referenced - declared
        if missing:
            raise ValueError(f"report references undeclared artifacts: {sorted(missing)}")
        return self


class ReportDelivery(BaseModel):
    model_config = ConfigDict(extra="forbid")

    report_path: Path
    manifest_path: Path
    report_sha256: str


def _metric(value: int | float | None, artifact_id: str) -> str:
    rendered = "NA" if value is None else str(value)
    return f"{rendered} <!--METRIC artifact:{artifact_id}-->"


def _environment() -> Environment:
    environment = Environment(
        loader=FileSystemLoader(TEMPLATE_ROOT),
        undefined=StrictUndefined,
        autoescape=select_autoescape(default_for_string=False),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    environment.globals["metric"] = _metric
    return environment


def _scan_claims(text: str) -> list[str]:
    return [phrase for phrase in PROHIBITED_CLAIMS if phrase in text]


def find_untraced_metrics(report: str) -> list[str]:
    """Return malformed metric markers lacking an artifact identifier."""
    markers = re.findall(r"<!--METRIC\s+([^>]*)-->", report)
    return [marker for marker in markers if not re.fullmatch(r"artifact:[A-Za-z0-9._-]+", marker.strip())]


def render_report(context: dict[str, Any], *, release: bool = False) -> str:
    """Validate and render the standard Chinese Markdown report."""
    validated = ReportContext.model_validate(context)
    if release:
        prohibited = _scan_claims(validated.conclusion)
        if prohibited:
            raise ProhibitedClaimError(
                "release conclusion contains prohibited computation-only claim(s): "
                + ", ".join(prohibited)
            )
    rendered = _environment().get_template("target_fishing_report.md.j2").render(
        context=validated
    )
    malformed = find_untraced_metrics(rendered)
    if malformed:
        raise ArtifactIntegrityError(f"untraced report metrics: {malformed}")
    return rendered.rstrip() + "\n"


def _validate_artifacts(context: ReportContext, artifact_root: Path) -> None:
    for artifact in context.artifacts:
        path = artifact_root / artifact.path
        if not path.is_file():
            raise ArtifactIntegrityError(
                f"artifact {artifact.artifact_id} does not exist: {path}"
            )
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed != artifact.sha256:
            raise ArtifactIntegrityError(
                f"artifact {artifact.artifact_id} hash mismatch: "
                f"expected {artifact.sha256}, observed {observed}"
            )


def write_report_delivery(
    context: dict[str, Any], output_dir: Path, *, artifact_root: Path
) -> ReportDelivery:
    """Validate all evidence artifacts and atomically create a delivery shell."""
    validated = ReportContext.model_validate(context)
    _validate_artifacts(validated, artifact_root)
    report = render_report(context, release=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    for directory in (
        "tables",
        "evidence_cards",
        "structures",
        "pymol",
        "trajectories",
        "audit",
    ):
        (output_dir / directory).mkdir(exist_ok=True)

    report_path = output_dir / "report.md"
    report_sha256 = write_bytes_atomic(report_path, report.encode("utf-8"))
    manifest = {
        "schema_version": "1.0",
        "project_id": validated.project_id,
        "release_status": "passed",
        "report": {"path": "report.md", "sha256": report_sha256},
        "template_sha256": content_sha256(
            (TEMPLATE_ROOT / "target_fishing_report.md.j2").read_bytes()
        ),
        "artifacts": [artifact.model_dump(mode="json") for artifact in validated.artifacts],
    }
    manifest_path = output_dir / "report_manifest.json"
    write_bytes_atomic(
        manifest_path,
        (json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode(),
    )
    return ReportDelivery(
        report_path=report_path,
        manifest_path=manifest_path,
        report_sha256=report_sha256,
    )


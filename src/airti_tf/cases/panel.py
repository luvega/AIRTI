"""Deterministic construction of a blinded engineering-panel target manifest."""

from __future__ import annotations

import hashlib
from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from airti_tf.cases.evaluation import load_case_definition
from airti_tf.manifest_io import read_jsonl, write_artifact, write_jsonl_atomic


class PanelSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: str
    panel_id: str
    panel_size: int = Field(gt=0)
    anchor_count: int = Field(gt=0)
    family_comparators_per_anchor: int = Field(gt=0)
    stratified_comparator_count: int = Field(gt=0)
    selection_seed: int
    strata: dict[str, int]
    family_comparators: dict[str, list[str]] = Field(default_factory=dict)
    comparison_proteins_are_true_negatives: bool

    @model_validator(mode="after")
    def validate_counts(self) -> "PanelSpec":
        expected = (
            self.anchor_count
            + self.anchor_count * self.family_comparators_per_anchor
            + self.stratified_comparator_count
        )
        if expected != self.panel_size:
            raise ValueError("panel count components do not equal panel_size")
        if sum(self.strata.values()) != self.stratified_comparator_count:
            raise ValueError("stratum quotas do not equal stratified comparator count")
        for anchor, comparators in self.family_comparators.items():
            if len(comparators) != self.family_comparators_per_anchor:
                raise ValueError(
                    f"frozen comparator count does not match quota for {anchor}"
                )
            if len(set(comparators)) != len(comparators):
                raise ValueError(f"duplicate frozen family comparators for {anchor}")
        return self


class CasePanelSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    panel_id: str
    target_count: int = Field(gt=0)
    row_count: int = Field(gt=0)
    manifest_path: Path
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    audit_path: Path
    audit_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


def _stable_key(target_id: str, *, seed: int) -> tuple[str, str]:
    material = f"{seed}:{target_id}".encode()
    return hashlib.sha256(material).hexdigest(), target_id


def _stratum(family: str) -> str:
    normalized = family.lower()
    if any(
        token in normalized
        for token in ("chaperone", "scaffold", "heat shock", "cytoskeleton")
    ):
        return "chaperone_or_scaffold"
    if any(token in normalized for token in ("transporter", "channel", "solute carrier")):
        return "transporter_or_channel"
    if any(token in normalized for token in ("receptor", "gpcr")):
        return "receptor"
    if any(
        token in normalized
        for token in (
            "enzyme",
            "kinase",
            "phosphatase",
            "protease",
            "peptidase",
            "cytochrome p450",
            "oxidase",
            "reductase",
            "transferase",
            "hydrolase",
        )
    ):
        return "enzyme"
    return "other"


def build_case_panel(
    *,
    case_path: Path,
    panel_spec_path: Path,
    target_manifest: Path,
    output_manifest: Path,
    audit_output: Path,
) -> CasePanelSummary:
    """Select exactly 64 ready targets without labeling comparators as negatives."""
    case = load_case_definition(case_path)
    raw_spec = yaml.safe_load(panel_spec_path.read_text(encoding="utf-8"))
    spec = PanelSpec.model_validate(raw_spec)
    anchors = [anchor.target_id for anchor in case.positive_anchors]
    if len(anchors) != spec.anchor_count:
        raise ValueError("case anchor count does not match panel specification")
    rows = read_jsonl(target_manifest)
    rows_by_target: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        target_id = str(row.get("target_id") or row.get("uniprot_id") or "")
        if target_id:
            rows_by_target.setdefault(target_id, []).append(row)
    status_by_target = {
        target_id: {
            str(
                row.get("status")
                or ("ready" if row.get("reviewed") is True else "unsupported")
            )
            for row in target_rows
        }
        for target_id, target_rows in rows_by_target.items()
    }
    ready_ids = {
        target_id
        for target_id, statuses in status_by_target.items()
        if statuses == {"ready"}
    }
    missing_anchors = sorted(set(anchors) - ready_ids)
    if missing_anchors:
        raise ValueError(f"panel anchors are not ready: {missing_anchors}")
    family_by_target = {
        target_id: str(
            target_rows[0].get("family")
            or target_rows[0].get("protein_family")
            or "unclassified"
        )
        for target_id, target_rows in rows_by_target.items()
    }
    selected = set(anchors)
    family_comparators: dict[str, list[str]] = {}
    for anchor in anchors:
        if anchor in spec.family_comparators:
            chosen = spec.family_comparators[anchor]
            unavailable = sorted(set(chosen) - ready_ids)
            overlapping = sorted(set(chosen) & selected)
            if unavailable:
                raise ValueError(
                    f"frozen family comparators are not ready for {anchor}: "
                    f"{unavailable}"
                )
            if overlapping:
                raise ValueError(
                    f"frozen family comparators overlap prior selections for "
                    f"{anchor}: {overlapping}"
                )
        else:
            family = family_by_target[anchor]
            pool = sorted(
                (
                    target_id
                    for target_id in ready_ids
                    if target_id not in selected
                    and family_by_target[target_id] == family
                ),
                key=lambda target_id: _stable_key(
                    target_id, seed=spec.selection_seed
                ),
            )
            chosen = pool[: spec.family_comparators_per_anchor]
            if len(chosen) != spec.family_comparators_per_anchor:
                raise ValueError(f"insufficient family comparators for {anchor}")
        family_comparators[anchor] = chosen
        selected.update(chosen)

    stratified: dict[str, list[str]] = {}
    for stratum, quota in spec.strata.items():
        pool = sorted(
            (
                target_id
                for target_id in ready_ids
                if target_id not in selected
                and _stratum(family_by_target[target_id]) == stratum
            ),
            key=lambda target_id: _stable_key(
                target_id, seed=spec.selection_seed
            ),
        )
        chosen = pool[:quota]
        if len(chosen) != quota:
            raise ValueError(f"insufficient ready targets for stratum {stratum}")
        stratified[stratum] = chosen
        selected.update(chosen)
    if len(selected) != spec.panel_size:
        raise RuntimeError("internal panel selection count mismatch")

    output_rows = [
        row
        for target_id in sorted(selected)
        for row in rows_by_target[target_id]
    ]
    manifest_sha256 = write_jsonl_atomic(output_manifest, output_rows)
    audit = {
        "schema_version": "1.0",
        "panel_id": spec.panel_id,
        "selection_seed": spec.selection_seed,
        "target_count": len(selected),
        "anchors": anchors,
        "family_comparators": family_comparators,
        "family_comparator_selection": (
            "frozen_ids" if spec.family_comparators else "exact_family"
        ),
        "strata": stratified,
        "comparison_proteins_are_true_negatives": (
            spec.comparison_proteins_are_true_negatives
        ),
        "source_manifest_sha256": hashlib.sha256(
            target_manifest.read_bytes()
        ).hexdigest(),
    }
    written = write_artifact(audit_output, audit)
    return CasePanelSummary(
        panel_id=spec.panel_id,
        target_count=len(selected),
        row_count=len(output_rows),
        manifest_path=output_manifest,
        manifest_sha256=manifest_sha256,
        audit_path=audit_output,
        audit_sha256=written.sha256,
    )

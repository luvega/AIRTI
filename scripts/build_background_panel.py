#!/usr/bin/env python3
"""Build or validate the immutable ChEMBL 37 background probe panel."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any, cast
from urllib.error import URLError
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen

import rdkit

from airti_tf.manifest_io import content_sha256, write_artifact
from airti_tf.screening.background import (
    BackgroundCandidate,
    build_chembl_candidate,
    select_background_panel,
    validate_panel_tsv,
    write_panel_tsv,
)

STATUS_URL = "https://www.ebi.ac.uk/chembl/api/data/status.json"
MOLECULE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
EXPECTED_RELEASE = "ChEMBL_37"
LICENSE = "Creative Commons Attribution-ShareAlike 3.0 Unported"


def fetch_json(url: str, *, attempts: int = 3) -> dict[str, Any]:
    request = Request(url, headers={"User-Agent": "AIRTI-Target-Fishing/0.1.0"})
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            with urlopen(request, timeout=120) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise ValueError(f"JSON root is not an object: {url}")
            return cast(dict[str, Any], payload)
        except (URLError, TimeoutError) as error:
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(f"failed to fetch {url}") from last_error


def fetch_chembl37_source() -> tuple[dict[str, Any], list[dict[str, str]]]:
    status = fetch_json(STATUS_URL)
    if status.get("chembl_db_version") != EXPECTED_RELEASE:
        raise RuntimeError(
            f"expected {EXPECTED_RELEASE}, API reports {status.get('chembl_db_version')}"
        )
    query = urlencode(
        {
            "max_phase": 4,
            "molecule_type": "Small molecule",
            "molecule_structures__isnull": "false",
            "limit": 1000,
        }
    )
    next_url: str | None = f"{MOLECULE_URL}?{query}"
    source: list[dict[str, str]] = []
    while next_url:
        page = fetch_json(next_url)
        molecules = page.get("molecules")
        if not isinstance(molecules, list):
            raise RuntimeError("ChEMBL molecule response lacks molecules list")
        for raw in molecules:
            if not isinstance(raw, dict):
                continue
            structures = raw.get("molecule_structures")
            if not isinstance(structures, dict):
                continue
            source_id = raw.get("molecule_chembl_id")
            smiles = structures.get("canonical_smiles")
            if isinstance(source_id, str) and isinstance(smiles, str):
                source.append({"source_id": source_id, "canonical_smiles": smiles})
        page_meta = page.get("page_meta")
        if not isinstance(page_meta, dict):
            raise RuntimeError("ChEMBL molecule response lacks page_meta")
        relative_next = page_meta.get("next")
        next_url = (
            urljoin("https://www.ebi.ac.uk", relative_next)
            if isinstance(relative_next, str)
            else None
        )
    return status, sorted(source, key=lambda row: row["source_id"])


def build_candidates(source: list[dict[str, str]]) -> list[BackgroundCandidate]:
    unique: dict[str, BackgroundCandidate] = {}
    for row in source:
        candidate = build_chembl_candidate(row["source_id"], row["canonical_smiles"])
        if candidate is None:
            continue
        previous = unique.get(candidate.canonical_smiles)
        if previous is None or candidate.source_id < previous.source_id:
            unique[candidate.canonical_smiles] = candidate
    return sorted(unique.values(), key=lambda item: item.source_id)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", type=Path)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument(
        "--output", type=Path, default=Path("data/reference/background_probes_v1.smi")
    )
    parser.add_argument(
        "--source-cache",
        type=Path,
        default=Path("data/reference/background_probes_v1.source.json"),
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_args()
    if arguments.check_only is not None:
        summary = validate_panel_tsv(arguments.check_only)
        print(
            f"valid background panel: {summary['probe_count']} probes, "
            f"{summary['strata']} strata, {summary['per_stratum']} per stratum"
        )
        return 0

    source_cache: Path = arguments.source_cache
    if arguments.fetch:
        status, source = fetch_chembl37_source()
        source_artifact = write_artifact(
            source_cache,
            {
                "schema_version": "1.0",
                "source": "ChEMBL approved small molecules",
                "source_url": MOLECULE_URL,
                "status_url": STATUS_URL,
                "release": status["chembl_db_version"],
                "release_date": status["chembl_release_date"],
                "license": LICENSE,
                "records": source,
            },
        )
        source_sha256 = source_artifact.sha256
    else:
        if not source_cache.is_file():
            raise SystemExit("source cache missing; run with --fetch")
        cached = json.loads(source_cache.read_text(encoding="utf-8"))
        if cached.get("release") != EXPECTED_RELEASE:
            raise SystemExit(f"source cache is not {EXPECTED_RELEASE}")
        status = {
            "chembl_db_version": cached["release"],
            "chembl_release_date": cached["release_date"],
        }
        source = cached["records"]
        source_sha256 = content_sha256(source_cache.read_bytes())

    candidates = build_candidates(source)
    panel = select_background_panel(
        candidates,
        panel_size=100,
        strata=10,
        seed=20_260_812,
    )
    output: Path = arguments.output
    panel_sha256 = write_panel_tsv(output, panel)
    validate_panel_tsv(output)
    metadata_path = output.with_suffix(".meta.json")
    write_artifact(
        metadata_path,
        {
            "schema_version": "1.0",
            "panel_version": "v1",
            "probe_count": 100,
            "mw_strata": 10,
            "probes_per_stratum": 10,
            "selection": "ECFP4 greedy MaxMin within molecular-weight quantile strata",
            "seed": 20_260_812,
            "source_release": status["chembl_db_version"],
            "source_release_date": status["chembl_release_date"],
            "source_url": MOLECULE_URL,
            "source_license": LICENSE,
            "source_record_count": len(source),
            "eligible_unique_candidate_count": len(candidates),
            "source_sha256": source_sha256,
            "panel_sha256": panel_sha256,
            "rdkit_version": rdkit.__version__,
            "builder_sha256": content_sha256(Path(__file__).read_bytes()),
        },
    )
    print(f"wrote {len(panel)} probes to {output} ({panel_sha256})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

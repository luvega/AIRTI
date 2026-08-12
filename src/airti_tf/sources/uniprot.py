"""UniProt human reference-proteome acquisition and normalization."""

from __future__ import annotations

import csv
import gzip
import hashlib
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, Self, cast
from urllib.parse import urlencode
from urllib.request import urlopen

from pydantic import BaseModel, ConfigDict, Field

from airti_tf.manifest_io import write_artifact, write_bytes_atomic, write_jsonl_atomic

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"


class UniProtRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: str = "1.0"
    uniprot_id: str
    gene_primary: str | None = None
    taxonomy_id: int
    sequence: str = Field(pattern=r"^[A-Z]+$")
    sequence_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    reviewed: bool
    release: str
    isoform_aliases: list[str] = Field(default_factory=list)


@dataclass(frozen=True)
class ProteomeSnapshot:
    manifest_path: Path
    summary_path: Path
    manifest_sha256: str


class BinaryResponse(Protocol):
    headers: Mapping[str, str]

    def read(self) -> bytes: ...

    def __enter__(self) -> Self: ...

    def __exit__(self, *args: object) -> object: ...


def build_uniprot_url(*, proteome: str) -> str:
    """Build an unfiltered reference-proteome stream URL."""
    parameters = {
        "compressed": "true",
        "format": "tsv",
        "fields": "accession,reviewed,gene_primary,organism_id,sequence",
        "query": f"proteome:{proteome}",
    }
    return f"{UNIPROT_STREAM_URL}?{urlencode(parameters)}"


def parse_uniprot_tsv(path: Path, *, release: str) -> list[UniProtRecord]:
    """Parse human canonical records while retaining isoforms as aliases."""
    candidates: dict[str, UniProtRecord] = {}
    aliases: dict[str, set[str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Entry", "Reviewed", "Organism (ID)", "Sequence"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"UniProt TSV lacks required columns: {sorted(required)}")
        for row in reader:
            accession = row["Entry"].strip()
            if not accession or int(row["Organism (ID)"]) != 9606:
                continue
            if "-" in accession:
                base_accession = accession.split("-", maxsplit=1)[0]
                aliases.setdefault(base_accession, set()).add(accession)
                continue
            sequence = "".join(row["Sequence"].split()).upper()
            if not sequence.isalpha() or not sequence.isascii():
                raise ValueError(f"invalid sequence alphabet for {accession}")
            reviewed = row["Reviewed"].strip().lower() == "reviewed"
            record = UniProtRecord(
                uniprot_id=accession,
                gene_primary=(row.get("Gene Names (primary)") or "").strip() or None,
                taxonomy_id=9606,
                sequence=sequence,
                sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                reviewed=reviewed,
                release=release,
            )
            previous = candidates.get(accession)
            if previous is None or (record.reviewed and not previous.reviewed):
                candidates[accession] = record

    records: list[UniProtRecord] = []
    for accession in sorted(candidates):
        records.append(
            candidates[accession].model_copy(
                update={"isoform_aliases": sorted(aliases.get(accession, set()))}
            )
        )
    return records


def write_proteome_snapshot(
    records: list[UniProtRecord],
    *,
    output_dir: Path,
    release: str,
    source_sha256: str,
    request_url: str,
    etag: str | None,
) -> ProteomeSnapshot:
    """Write the canonical JSONL manifest and its audit summary."""
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "human_canonical_proteome.jsonl"
    manifest_sha256 = write_jsonl_atomic(
        manifest_path,
        [record.model_dump(mode="json") for record in records],
    )
    summary = {
        "schema_version": "1.0",
        "release": release,
        "request_url": request_url,
        "etag": etag,
        "source_sha256": source_sha256,
        "manifest_sha256": manifest_sha256,
        "record_count": len(records),
        "reviewed_count": sum(record.reviewed for record in records),
        "unreviewed_count": sum(not record.reviewed for record in records),
        "taxonomy_id": 9606,
    }
    summary_path = output_dir / "human_canonical_proteome.summary.json"
    write_artifact(summary_path, summary)
    return ProteomeSnapshot(
        manifest_path=manifest_path,
        summary_path=summary_path,
        manifest_sha256=manifest_sha256,
    )


def _default_opener(url: str) -> BinaryResponse:
    return cast(BinaryResponse, urlopen(url, timeout=120))


def fetch_uniprot_snapshot(
    *,
    proteome: str,
    release: str,
    output_dir: Path,
    opener: Callable[[str], BinaryResponse] = _default_opener,
) -> ProteomeSnapshot:
    """Download a release-pinned raw snapshot and build its canonical manifest."""
    request_url = build_uniprot_url(proteome=proteome)
    with opener(request_url) as response:
        compressed = response.read()
        response_release = response.headers.get("X-UniProt-Release")
        etag = response.headers.get("ETag")
    if response_release != release:
        raise ValueError(
            f"requested UniProt release {release}, received {response_release or 'unknown'}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    compressed_path = output_dir / f"uniprot_human_{release}.tsv.gz"
    source_sha256 = write_bytes_atomic(compressed_path, compressed)
    raw_path = output_dir / f"uniprot_human_{release}.tsv"
    try:
        raw_tsv = gzip.decompress(compressed)
    except gzip.BadGzipFile as error:
        raise ValueError("UniProt response is not valid gzip data") from error
    write_bytes_atomic(raw_path, raw_tsv)
    records = parse_uniprot_tsv(raw_path, release=release)
    return write_proteome_snapshot(
        records,
        output_dir=output_dir,
        release=release,
        source_sha256=source_sha256,
        request_url=request_url,
        etag=etag,
    )

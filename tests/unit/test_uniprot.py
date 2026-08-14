import gzip
import io
import json
from pathlib import Path

from airti_tf.sources.uniprot import (
    build_uniprot_url,
    fetch_uniprot_snapshot,
    parse_uniprot_tsv,
    write_proteome_snapshot,
)


def test_builds_canonical_manifest() -> None:
    fixture = Path("tests/fixtures/uniprot_human_sample.tsv")

    records = parse_uniprot_tsv(fixture, release="2026_03")

    assert [record.uniprot_id for record in records] == ["P00533", "P04637"]
    assert all(record.taxonomy_id == 9606 for record in records)
    assert all(record.sequence_sha256 for record in records)
    assert records[0].reviewed is True
    assert records[0].isoform_aliases == ["P00533-2"]


def test_production_query_is_reviewed_human_reference_proteome() -> None:
    url = build_uniprot_url(proteome="UP000005640")

    assert "proteome%3AUP000005640" in url
    assert "reviewed%3Atrue" in url
    assert "organism_id%3A9606" in url
    assert "ft_transmem" in url
    assert "cc_subcellular_location" in url


def test_transmembrane_annotations_are_counted(tmp_path: Path) -> None:
    snapshot = tmp_path / "uniprot.tsv"
    snapshot.write_text(
        "Entry\tReviewed\tGene Names (primary)\tOrganism (ID)\tSequence\t"
        "Transmembrane\n"
        "P08183\treviewed\tABCB1\t9606\tMPEPTIDE\t"
        'TRANSMEM 10..30; /note="Helical"; TRANSMEM 50..70; /note="Helical"\n',
        encoding="utf-8",
    )

    records = parse_uniprot_tsv(snapshot, release="2026_03")

    assert records[0].transmembrane_segments == 2


def test_membrane_association_is_inferred_from_subcellular_annotation(
    tmp_path: Path,
) -> None:
    snapshot = tmp_path / "uniprot.tsv"
    snapshot.write_text(
        "Entry\tReviewed\tGene Names (primary)\tOrganism (ID)\tSequence\t"
        "Subcellular location [CC]\tProtein families\n"
        "P20815\treviewed\tCYP3A5\t9606\tMPEPTIDE\t"
        "SUBCELLULAR LOCATION: Endoplasmic reticulum membrane; Peripheral "
        "membrane protein.\tCytochrome P450 family\n"
        "P07900\treviewed\tHSP90AA1\t9606\tMPEPTIDE\t"
        "SUBCELLULAR LOCATION: Cell membrane; Peripheral membrane protein.\t"
        "Heat shock protein 90 family\n",
        encoding="utf-8",
    )

    records = parse_uniprot_tsv(snapshot, release="2026_03")

    by_id = {record.uniprot_id: record for record in records}
    assert by_id["P20815"].transmembrane_segments == 0
    assert by_id["P20815"].membrane_associated is True
    assert by_id["P07900"].membrane_associated is False


def test_unreviewed_only_entries_are_excluded_from_production_snapshot() -> None:
    records = parse_uniprot_tsv(
        Path("tests/fixtures/uniprot_human_sample.tsv"), release="2026_03"
    )

    assert all(record.reviewed for record in records)
    assert "QUNREV" not in {record.uniprot_id for record in records}


def test_snapshot_is_versioned_and_summarized(tmp_path: Path) -> None:
    records = parse_uniprot_tsv(
        Path("tests/fixtures/uniprot_human_sample.tsv"), release="2026_03"
    )

    snapshot = write_proteome_snapshot(
        records,
        output_dir=tmp_path,
        release="2026_03",
        source_sha256="f" * 64,
        request_url="https://rest.uniprot.org/example",
        etag='"fixture"',
    )

    assert snapshot.manifest_path.name == "human_canonical_proteome.jsonl"
    summary = json.loads(snapshot.summary_path.read_text(encoding="utf-8"))
    assert summary["record_count"] == 2
    assert summary["reviewed_count"] == 2
    assert summary["unreviewed_count"] == 0
    assert summary["release"] == "2026_03"
    assert summary["source_sha256"] == "f" * 64


def test_fetch_records_release_etag_and_raw_hash(tmp_path: Path) -> None:
    raw_tsv = Path("tests/fixtures/uniprot_human_sample.tsv").read_bytes()
    compressed_tsv = gzip.compress(raw_tsv, mtime=0)

    class Response(io.BytesIO):
        headers = {"X-UniProt-Release": "2026_03", "ETag": '"fixture-etag"'}

    def opener(request_url: str) -> Response:
        assert "UP000005640" in request_url
        return Response(compressed_tsv)

    snapshot = fetch_uniprot_snapshot(
        proteome="UP000005640",
        release="2026_03",
        output_dir=tmp_path,
        opener=opener,
    )

    summary = json.loads(snapshot.summary_path.read_text(encoding="utf-8"))
    assert summary["etag"] == '"fixture-etag"'
    assert len(summary["source_sha256"]) == 64
    assert (tmp_path / "uniprot_human_2026_03.tsv.gz").read_bytes() == compressed_tsv

from pathlib import Path
import shutil

import pytest

from airti_tf.refinement.boltz2 import (
    BoltzJob,
    BoltzExecutionResult,
    BoltzSeedResult,
    InsufficientBoltzSeedsError,
    MissingMSAError,
    build_boltz_command,
    build_boltz_yaml,
    classify_boltz_failure,
    msa_cache_path,
    parse_boltz_output,
    assess_boltz_structure,
    run_boltz_seed,
    summarize_boltz_seeds,
    write_boltz_quality,
)


@pytest.fixture
def msa(tmp_path: Path) -> Path:
    path = tmp_path / "P00533.a3m"
    path.write_text(">query\nMPEPTIDE\n", encoding="utf-8")
    return path


@pytest.fixture
def job(msa: Path, tmp_path: Path) -> BoltzJob:
    return BoltzJob(
        job_id="job1",
        target_id="P00533",
        sequence="MPEPTIDE",
        sequence_sha256="a" * 64,
        msa_path=msa,
        msa_database_version="af3-db-2022-2023",
        ligand_state_id="L1:S1",
        ligand_smiles="CCO",
        ligand_atom_count=3,
        pocket_residues=[718, 719],
        input_yaml=tmp_path / "job1.yaml",
        output_dir=tmp_path / "seed11",
        cache_path=tmp_path / "boltz-cache",
    )


def test_boltz_yaml_contains_affinity_and_pocket_constraint(job: BoltzJob) -> None:
    payload = build_boltz_yaml(job)

    assert payload["properties"][0]["affinity"]["binder"] == "B"
    assert payload["constraints"][0]["pocket"]["binder"] == "B"
    assert payload["constraints"][0]["pocket"]["contacts"] == [["A", 718], ["A", 719]]
    assert payload["sequences"][0]["protein"]["msa"].endswith(".a3m")


def test_missing_cached_msa_blocks_production(job: BoltzJob) -> None:
    missing = job.model_copy(update={"msa_path": Path("missing.a3m")})

    with pytest.raises(MissingMSAError):
        build_boltz_yaml(missing, profile="production")


def test_msa_cache_identity_includes_sequence_and_database() -> None:
    path = msa_cache_path(
        Path("/cache/msa"),
        uniprot_id="P00533",
        sequence_sha256="a" * 64,
        database_version="af3-db-2022-2023",
    )

    assert path.name == f"P00533.{'a' * 12}.af3-db-2022-2023.a3m"


def test_command_uses_three_samples_seed_and_potentials(job: BoltzJob) -> None:
    command = build_boltz_command(job, seed=29)

    assert command[:3] == ["boltz", "predict", str(job.input_yaml)]
    assert command[command.index("--diffusion_samples") + 1] == "3"
    assert command[command.index("--diffusion_samples_affinity") + 1] == "3"
    assert command[command.index("--seed") + 1] == "29"
    assert command[command.index("--cache") + 1] == str(job.cache_path)
    assert "--use_potentials" in command


def test_parser_reads_official_confidence_and_affinity_outputs() -> None:
    result = parse_boltz_output(
        Path("tests/fixtures/boltz2_output"), input_stem="job1", seed=11
    )

    assert result.status == "succeeded"
    assert result.confidence_score == pytest.approx(0.81)
    assert result.ligand_iptm == pytest.approx(0.74)
    assert result.affinity_probability == pytest.approx(0.86)
    assert result.affinity_pred_value == pytest.approx(-1.3)
    assert result.pocket_constraint_fraction == pytest.approx(0.92)
    assert result.severe_clash is False


def successful(seed: int, confidence: float, affinity: float) -> BoltzSeedResult:
    return BoltzSeedResult(
        seed=seed,
        status="succeeded",
        confidence_score=confidence,
        ligand_iptm=confidence - 0.05,
        affinity_probability=affinity,
        affinity_pred_value=-1.0,
        pocket_constraint_fraction=0.9,
        severe_clash=False,
    )


def test_two_of_three_seeds_form_median_consensus() -> None:
    summary = summarize_boltz_seeds(
        [
            successful(11, 0.8, 0.9),
            successful(29, 0.7, 0.7),
            BoltzSeedResult(seed=47, status="failed", error_code="cuda_oom"),
        ]
    )

    assert summary.seed_success_count == 2
    assert summary.confidence_median == pytest.approx(0.75)
    assert summary.affinity_probability_median == pytest.approx(0.8)


def test_clashing_seed_does_not_count_as_success() -> None:
    clashing = successful(47, 0.9, 0.9).model_copy(update={"severe_clash": True})

    with pytest.raises(InsufficientBoltzSeedsError):
        summarize_boltz_seeds(
            [successful(11, 0.8, 0.8), clashing, BoltzSeedResult(seed=29, status="failed", error_code="nan_output")]
        )


@pytest.mark.parametrize(
    ("stderr", "expected"),
    [
        ("CUDA out of memory", "cuda_oom"),
        ("invalid yaml", "invalid_yaml"),
        ("ligand has 129 atoms", "ligand_too_large"),
        ("nan in output", "nan_output"),
        ("pocket constraint violation", "constraint_violation"),
    ],
)
def test_failure_classification(stderr: str, expected: str) -> None:
    assert classify_boltz_failure(stderr) == expected


def test_numba_cache_trace_is_not_misclassified_as_missing_msa() -> None:
    stderr = (
        "File '/opt/conda/lib/python3.11/site-packages/boltz/data/feature/"
        "featurizer.py', in _prepare_msa_arrays_inner\n"
        "RuntimeError: cannot cache function: no locator available"
    )

    assert classify_boltz_failure(stderr) == "runtime_environment"


def test_structural_qc_measures_pocket_contacts_and_clashes(tmp_path: Path) -> None:
    structure = tmp_path / "prediction.cif"
    structure.write_text(
        """data_model
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
ATOM 1 C CA ALA 1 1 A 0.0 0.0 0.0 A
ATOM 2 C CA GLY 2 2 A 20.0 0.0 0.0 A
HETATM 3 C C1 LIG1 . . B 3.0 0.0 0.0 B
#
""",
        encoding="utf-8",
    )

    quality = assess_boltz_structure(structure, pocket_residues=[1, 2])

    assert quality.pocket_constraint_fraction == pytest.approx(0.5)
    assert quality.severe_clash is False
    quality_path = write_boltz_quality(structure, pocket_residues=[1, 2])
    assert quality_path.name == "airti_quality.json"
    assert quality_path.is_file()


def test_boltz_runner_locates_official_result_directory(job: BoltzJob) -> None:
    def executor(command: list[str], timeout_seconds: int) -> BoltzExecutionResult:
        assert timeout_seconds == 7200
        output_dir = Path(command[command.index("--out_dir") + 1])
        result_root = output_dir / f"boltz_results_{job.input_yaml.stem}"
        shutil.copytree("tests/fixtures/boltz2_output/predictions", result_root / "predictions")
        source = result_root / "predictions/job1"
        destination = result_root / "predictions" / job.input_yaml.stem
        source.rename(destination)
        for path in destination.iterdir():
            path.rename(path.with_name(path.name.replace("job1", job.input_yaml.stem)))
        return BoltzExecutionResult(
            return_code=0,
            stdout="Number of failed examples: 0",
            stderr="",
            timed_out=False,
        )

    result = run_boltz_seed(job, seed=11, executor=executor)

    assert result.status == "succeeded"
    assert result.confidence_score == pytest.approx(0.81)
    assert result.structure_path is not None
    assert result.structure_path.is_file()
    assert (job.output_dir / "boltz.stdout.log").read_text() == (
        "Number of failed examples: 0"
    )
    assert (job.output_dir / "boltz.stderr.log").read_text() == ""
    assert (job.output_dir / "boltz.execution.json").is_file()

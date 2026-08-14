from urllib.error import URLError

from airti_tf.pockets.fpocket import PocketCandidate, PocketQC
from airti_tf.targets.production import (
    _clean_structure,
    _patch_heme_cif_charges,
    _read_url,
    _select_dockable_pockets,
)


def _atom_line(*, serial: int, altloc: str, occupancy: float) -> str:
    return (
        f"ATOM  {serial:5d}  CA {altloc}ALA A   1      "
        f"  1.000   2.000   3.000{occupancy:6.2f} 20.00           C  \n"
    )


def test_structure_cleanup_keeps_one_deterministic_altloc() -> None:
    raw = (
        "HEADER    TEST\n"
        + _atom_line(serial=1, altloc="A", occupancy=0.60)
        + _atom_line(serial=2, altloc="B", occupancy=0.40)
        + "END\n"
    ).encode()

    cleaned, retained_hem = _clean_structure(
        raw, chain_id="A", keep_hem=True
    )

    text = cleaned.decode()
    assert retained_hem is False
    assert text.count("ATOM  ") == 1
    assert " AALA" not in text
    assert " BALA" not in text


def test_https_transport_failure_is_retried_with_bounded_backoff(
    monkeypatch,
) -> None:
    attempts = 0
    sleeps: list[float] = []

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *_args: object) -> None:
            return None

        @staticmethod
        def read() -> bytes:
            return b"ok"

    def opener(_url: str, *, timeout: int):
        nonlocal attempts
        assert timeout == 120
        attempts += 1
        if attempts == 1:
            raise URLError("transient TLS EOF")
        return Response()

    monkeypatch.setattr("airti_tf.targets.production.urlopen", opener)
    monkeypatch.setattr(
        "airti_tf.targets.production.time.sleep", sleeps.append
    )

    assert _read_url("https://example.test/structure") == b"ok"
    assert attempts == 2
    assert sleeps == [0.5]


def test_heme_cif_patch_is_narrow_and_auditable() -> None:
    raw = """data_HEM
HEM FE FE FE ? 0.0 0.0 0.0
HEM O2A O2A O -1 1.0 1.0 1.0
HEM O2D O2D O -1 2.0 2.0 2.0
HEM CHA CHA C 0 3.0 3.0 3.0
HEM O2A H2A SING N N 77 N N
HEM FE NA SING N N 79 Y N
"""

    patched = _patch_heme_cif_charges(raw)

    assert "HEM FE FE FE 0 0.0 0.0 0.0" in patched
    assert "HEM O2A O2A O 0 1.0 1.0 1.0" in patched
    assert "HEM O2D O2D O 0 2.0 2.0 2.0" in patched
    assert "HEM CHA CHA C 0 3.0 3.0 3.0" in patched
    assert "HEM O2A H2A SING N N 77 N N" in patched
    assert "HEM FE NA SING N N 79 Y N" in patched


def test_heme_cif_patch_refuses_an_unexpected_component() -> None:
    raw = """data_HEM
HEM FE FE FE ? 0.0 0.0 0.0
HEM O2A O2A O -1 1.0 1.0 1.0
"""

    try:
        _patch_heme_cif_charges(raw)
    except ValueError as error:
        assert "O2D" in str(error)
    else:
        raise AssertionError("missing HEM atom must fail closed")


def _pocket_atom(x: float, y: float, z: float) -> str:
    return (
        f"ATOM      1  C   STP A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
    )


def _qualified(rank: int) -> PocketQC:
    return PocketQC(
        pocket=PocketCandidate(
            pocket_id=f"P00000:hash:{rank}",
            target_id="P00000",
            rank=rank,
            volume_a3=200,
            druggability=0.8,
            fpocket_score=10,
            residue_count=10,
        ),
        status="ready",
        score=0.8,
    )


def test_oversized_top_pocket_falls_through_to_next_candidate(tmp_path) -> None:
    pockets = tmp_path / "pockets"
    pockets.mkdir()
    (pockets / "pocket1_atm.pdb").write_text(
        _pocket_atom(0, 0, 0) + _pocket_atom(40, 0, 0),
        encoding="utf-8",
    )
    (pockets / "pocket2_atm.pdb").write_text(
        _pocket_atom(0, 0, 0) + _pocket_atom(5, 5, 5),
        encoding="utf-8",
    )

    selected = _select_dockable_pockets(
        [_qualified(1), _qualified(2)],
        pocket_output=tmp_path,
        max_pockets=1,
    )

    assert len(selected) == 1
    assert selected[0][0].pocket.rank == 2

import json
from pathlib import Path

import pytest

from airti_tf.reporting.render import (
    ProhibitedClaimError,
    find_untraced_metrics,
    render_report,
)


@pytest.fixture
def report_context() -> dict[str, object]:
    return json.loads(Path("tests/fixtures/report_context.json").read_text())


def test_every_report_metric_has_artifact_provenance(
    report_context: dict[str, object],
) -> None:
    report = render_report(report_context)

    assert find_untraced_metrics(report) == []
    assert "artifact:coverage-table" in report
    assert "artifact:rank-table" in report


@pytest.mark.parametrize("phrase", ["确认靶点", "已证实直接结合", "真实结合概率"])
def test_prohibited_claims_block_release(
    report_context: dict[str, object], phrase: str
) -> None:
    report_context["conclusion"] = f"本研究{phrase} EGFR"

    with pytest.raises(ProhibitedClaimError):
        render_report(report_context, release=True)


def test_report_has_required_scientific_boundaries(
    report_context: dict[str, object],
) -> None:
    report = render_report(report_context, release=True)

    assert "候选靶点" in report
    assert "unsupported" in report
    assert "湿实验建议" in report
    assert "局限性" in report


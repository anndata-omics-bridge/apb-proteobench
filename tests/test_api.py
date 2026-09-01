"""ProteoBench annotation, APB persistence, protocol, and CLI behavior."""

from __future__ import annotations

from pathlib import Path

import pytest
from apb2.result_facade import JsonValue, read_parsed_levels, write_parsed_levels
from loguru import logger

from apb_proteobench.annotation import ProteoBenchAnnotationParser
from apb_proteobench.api import annotate_result, score_result
from apb_proteobench.calculation.contracts import QuantitativeLevelInput
from apb_proteobench.calculation.intermediate import IntermediateResult
from apb_proteobench.calculation.metrics import ProteoBenchScores
from apb_proteobench.cli import app
from apb_proteobench.configuration.schema import ModuleSettings
from apb_proteobench.workflow import (
    MixedSpeciesDiagnostics,
    ProteoBenchCompatibleScoring,
    analyze_level,
)
from conftest import module_settings, parsed_result, quantitative_input, write_module


def test_annotation_requires_exact_coverage_and_embeds_complete_configuration(
    tmp_path: Path,
) -> None:
    module = tmp_path / "module.toml"
    write_module(module, alias=True)
    parsed = parsed_result()

    result = ProteoBenchAnnotationParser.from_path(module).parse(parsed).annotate()

    level = result.parsed.levels["ion"]
    assert level.obs.frame.columns == ["Run", "raw_file", "sample_name", "condition"]
    assert level.obs.frame.get_column("sample_name").to_list() == ["A1", "A2", "B1", "B2"]
    record = result.parsed.metadata["annotation"]
    assert isinstance(record, dict)
    proteobench = record["proteobench"]
    assert isinstance(proteobench, dict)
    details = _object(proteobench["metadata"])
    source = _object(details["source"])
    configuration = _object(details["configuration"])
    general = _object(configuration["general"])
    assert source["sha256"]
    assert general["level"] == "ion"
    assert parsed.levels["ion"].obs.frame.columns == ["Run"]


def test_annotation_rejects_quantification_and_module_subsets(tmp_path: Path) -> None:
    module = tmp_path / "module.toml"
    write_module(module)
    parser = ProteoBenchAnnotationParser.from_path(module)
    missing_quant = parsed_result()
    missing_quant.levels["ion"].obs.frame = missing_quant.levels["ion"].obs.frame.head(3)
    layer = missing_quant.levels["ion"].layers["Intensity"]
    layer.values = layer.values.select(layer.values.columns[:-1])

    with pytest.raises(ValueError, match="samples absent from quantification"):
        parser.parse(missing_quant)

    missing_module = tmp_path / "missing.toml"
    text = module.read_text(encoding="utf-8")
    missing_module.write_text(text.rsplit("[[samples]]", maxsplit=1)[0], encoding="utf-8")
    with pytest.raises(ValueError, match="complete sample annotation required"):
        ProteoBenchAnnotationParser.from_path(missing_module).parse(parsed_result())


@pytest.mark.parametrize("suffix", [".h5ad", ".h5mu", ".parquet", ".duckdb"])
def test_annotation_scoring_and_roundtrip_through_every_apb_format(
    suffix: str,
    tmp_path: Path,
) -> None:
    source = tmp_path / f"source{suffix}"
    annotated = tmp_path / f"annotated{suffix}"
    scored = tmp_path / f"scored{suffix}"
    module = tmp_path / "module.toml"
    write_module(module)
    write_parsed_levels(parsed_result(), source)

    annotate_result(source, module, annotated)
    result = score_result(annotated, scored)
    restored = read_parsed_levels(scored)

    assert result.analysis.scores.nr_feature == 3
    assert restored.levels["ion"].varm["proteobench"].height == 6
    stored = restored.levels["ion"].metadata["proteobench"]
    assert isinstance(stored, dict)
    assert _object(stored["scores"])["nr_feature"] == 3
    assert _object(stored["column_roles"])["Proteins"] == "var:Protein_Ids"
    with pytest.raises(ValueError, match="already exists"):
        score_result(annotated, scored)


class _ObservedDiagnostics:
    def __init__(self, result: IntermediateResult) -> None:
        self.result = result
        self.called = False

    def diagnose(
        self,
        level: QuantitativeLevelInput,
        configuration: ModuleSettings,
        /,
    ) -> IntermediateResult:
        del level, configuration
        self.called = True
        return self.result

    def identity(self) -> dict[str, str]:
        return {"name": "test-diagnostics", "version": "1"}


class _ObservedScoring:
    def __init__(self, result: ProteoBenchScores) -> None:
        self.result = result
        self.called = False

    def score(
        self,
        diagnostics: IntermediateResult,
        configuration: ModuleSettings,
        /,
    ) -> ProteoBenchScores:
        del diagnostics, configuration
        self.called = True
        return self.result

    def identity(self) -> dict[str, str]:
        return {"name": "test-scoring", "version": "1"}


def test_protocol_implementations_are_substitutable() -> None:
    inputs = quantitative_input()
    configuration = module_settings()
    expected = analyze_level(
        inputs,
        configuration,
        MixedSpeciesDiagnostics(),
        ProteoBenchCompatibleScoring(),
    )
    diagnostics = _ObservedDiagnostics(expected.diagnostics)
    scoring = _ObservedScoring(expected.scores)

    result = analyze_level(inputs, configuration, diagnostics, scoring)

    assert diagnostics.called and scoring.called
    assert result.diagnostic_method["name"] == "test-diagnostics"
    assert result.scoring_method["name"] == "test-scoring"


def test_cli_annotate_score_and_verbose_summary(tmp_path: Path) -> None:
    source = tmp_path / "source.parquet"
    annotated = tmp_path / "annotated.parquet"
    scored = tmp_path / "scored.parquet"
    concise = tmp_path / "concise.parquet"
    module = tmp_path / "module.toml"
    write_module(module)
    write_parsed_levels(parsed_result(), source)
    messages: list[str] = []
    sink = logger.add(messages.append, format="{message}")
    try:
        assert (
            app(
                ["annotate", str(source), str(module), str(annotated)],
                exit_on_error=False,
                result_action="return_value",
            )
            == 0
        )
        assert (
            app(
                ["score", str(annotated), str(concise)],
                exit_on_error=False,
                result_action="return_value",
            )
            == 0
        )
        assert (
            app(
                ["score", str(annotated), str(scored), "--verbose"],
                exit_on_error=False,
                result_action="return_value",
            )
            == 0
        )
    finally:
        logger.remove(sink)

    rendered = "".join(messages)
    assert "ProteoBench score summary" in rendered
    assert "species=['YEAST', 'ECOLI', 'HUMAN']" in rendered
    assert "features total=6 included=3 excluded=3" in rendered
    assert "varm['proteobench']" in rendered
    assert "scored level=ion" in rendered


def _object(value: JsonValue) -> dict[str, JsonValue]:
    assert isinstance(value, dict)
    return value

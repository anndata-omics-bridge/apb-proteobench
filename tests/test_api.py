"""ProteoBench annotation, APB persistence, protocol, and CLI behavior."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import polars as pl
import pytest
from apb2.result_facade import (
    AnnotationTable,
    FeatureRelation,
    FinalLayerTable,
    JsonValue,
    ParsedLevels,
    read_parsed_levels,
    write_parsed_levels,
)
from loguru import logger

from apb_proteobench.annotation import ProteoBenchAnnotationParser
from apb_proteobench.api import annotate_result, benchmark_result, score_result
from apb_proteobench.calculation.contracts import QuantitativeLevelInput
from apb_proteobench.calculation.intermediate import IntermediateResult
from apb_proteobench.calculation.metrics import ProteoBenchScores
from apb_proteobench.cli import app
from apb_proteobench.configuration.schema import ModuleSettings
from apb_proteobench.integration import (
    ALL_ABUNDANCE_LAYERS,
    PRIMARY_LAYER,
    NamedAbundanceLayer,
)
from apb_proteobench.workflow import (
    MixedSpeciesDiagnostics,
    ProteoBenchCompatibleScoring,
    analyze_level,
)
from conftest import module_settings, parsed_result, quantitative_input, write_module


def test_public_api_uses_only_apb2_public_facades() -> None:
    source = (Path(__file__).parents[1] / "src/apb_proteobench/api.py").read_text(encoding="utf-8")

    assert "apb2.parserV2" not in source


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
    record = _object(result.parsed.metadata["proteobench"])
    details = _object(_object(record["provenance"])["annotation"])
    assert set(record) == {"provenance"}
    assert "metadata" not in details and "convention" not in details
    source = _object(details["source"])
    configuration = _object(details["configuration"])
    general = _object(configuration["general"])
    assert source["sha256"]
    assert general["level"] == "ion"
    assert parsed.levels["ion"].obs.frame.columns == ["Run"]


def test_annotation_accepts_module_sample_names_as_identifiers(tmp_path: Path) -> None:
    module = tmp_path / "module.toml"
    write_module(module)
    parsed = parsed_result()
    parsed.levels["ion"].obs.frame = parsed.levels["ion"].obs.frame.with_columns(
        pl.Series("Run", ["A1", "A2", "B1", "B2"])
    )

    result = ProteoBenchAnnotationParser.from_path(module).parse(parsed).annotate()

    level = result.parsed.levels["ion"]
    assert level.obs.frame.get_column("raw_file").to_list() == [
        "run_A1",
        "run_A2",
        "run_B1",
        "run_B2",
    ]
    assert level.obs.frame.get_column("sample_name").to_list() == ["A1", "A2", "B1", "B2"]


def test_annotation_preserves_root_annotation_tables_and_relations(tmp_path: Path) -> None:
    module = tmp_path / "module.toml"
    write_module(module, alias=True)
    parsed = parsed_result()
    parsed.annotation_tables["members"] = AnnotationTable(
        frame=pl.DataFrame({"member": ["P1"]}),
        key_columns=("member",),
    )
    parsed.feature_relations["membership"] = FeatureRelation(
        annotation_table="members",
        target_level="ion",
        coordinates=pl.DataFrame({"row": [0], "column": [0], "value": [1.0]}),
    )

    result = ProteoBenchAnnotationParser.from_path(module).parse(parsed).annotate().parsed

    assert result.annotation_tables["members"].frame.equals(
        parsed.annotation_tables["members"].frame
    )
    assert result.feature_relations["membership"].coordinates.equals(
        parsed.feature_relations["membership"].coordinates
    )


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
    baseline = read_parsed_levels(annotated)

    assert result.layers["Intensity"].analysis.scores.nr_feature == 3
    assert restored.levels["ion"].varm["proteobench:Intensity"].height == 6
    stored = restored.levels["ion"].metadata["proteobench"]
    assert isinstance(stored, dict)
    provenance = _object(_object(restored.metadata["proteobench"])["provenance"])
    assert _object(provenance["scoring"])["schema_version"] == "3"
    assert set(provenance) == {"annotation", "scoring"}
    assert (
        stored["annotation"]
        == _object(baseline.levels["ion"].metadata["proteobench"])["annotation"]
    )
    assert set(stored) == {"annotation", "scoring"}
    layer = _object(_object(stored["scoring"])["Intensity"])
    assert _object(layer["scores"])["nr_feature"] == 3
    assert _object(layer["column_roles"])["Proteins"] == "var:Protein_Ids"
    assert (
        provenance["annotation"]
        == _object(_object(baseline.metadata["proteobench"])["provenance"])["annotation"]
    )
    before, after = baseline.levels["ion"], restored.levels["ion"]
    assert before.uns == after.uns
    for first, second in [(before.obs.frame, after.obs.frame), (before.var.frame, after.var.frame)]:
        assert first.equals(second) and first.schema == second.schema
    for name, matrix in before.layers.items():
        assert matrix.values.equals(after.layers[name].values)
        assert matrix.values.schema == after.layers[name].values.schema
    for name, table in before.varm.items():
        assert table.equals(after.varm[name])
    with pytest.raises(ValueError, match="already exists"):
        score_result(annotated, scored)
    target = tmp_path / f"rescored{suffix}"
    with pytest.raises(ValueError, match="refusing to overwrite"):
        score_result(scored, target)
    assert not target.exists()


def test_benchmark_result_annotates_and_scores_in_one_call(tmp_path: Path) -> None:
    source = tmp_path / "source.parquet"
    target = tmp_path / "benchmarked.parquet"
    module = tmp_path / "module.toml"
    write_module(module)
    write_parsed_levels(parsed_result(), source)

    result = benchmark_result(source, module, target)

    restored = read_parsed_levels(target)
    assert result.layers["Intensity"].analysis.scores.nr_feature == 3
    assert "sample_name" in restored.levels["ion"].obs.frame
    assert "proteobench:Intensity" in restored.levels["ion"].varm


def _multi_layer_result() -> ParsedLevels:
    parsed = parsed_result()
    level = parsed.levels["ion"]
    values = level.layers["Intensity"].values
    level.layers["LFQ/Intensity"] = FinalLayerTable(
        layer_name="LFQ/Intensity",
        var_key_columns=("feature",),
        values=values.clone(),
    )
    level.layers["QValue"] = FinalLayerTable(
        layer_name="QValue",
        var_key_columns=("feature",),
        values=values.clone(),
    )
    level.uns["layer_roles"] = {"abundance": ["Intensity", "LFQ/Intensity"]}
    return parsed


def _write_annotated_result(
    folder: Path,
    parsed: ParsedLevels,
    /,
    *,
    suffix: str = ".parquet",
) -> Path:
    source = folder / f"source{suffix}"
    annotated = folder / f"annotated{suffix}"
    module = folder / "module.toml"
    write_module(module)
    write_parsed_levels(parsed, source)
    annotate_result(source, module, annotated)
    return annotated


def test_default_selection_scores_every_abundance_layer(tmp_path: Path) -> None:
    annotated = _write_annotated_result(tmp_path, _multi_layer_result())
    target = tmp_path / "scored.parquet"

    result = score_result(annotated, target)

    restored = read_parsed_levels(target)
    assert list(result.layers) == ["Intensity", "LFQ/Intensity"]
    assert set(restored.levels["ion"].varm) == {
        "proteobench:Intensity",
        "proteobench:LFQ/Intensity",
    }
    record = _object(restored.levels["ion"].metadata["proteobench"])
    provenance = _object(
        _object(_object(restored.metadata["proteobench"])["provenance"])["scoring"]
    )
    assert provenance["selection_mode"] == "all_abundance"
    assert "resolved_layers" not in provenance
    assert list(_object(record["scoring"])) == ["Intensity", "LFQ%2FIntensity"]


def test_primary_selection_scores_only_x_layer(tmp_path: Path) -> None:
    annotated = _write_annotated_result(tmp_path, _multi_layer_result())
    target = tmp_path / "scored.parquet"

    result = score_result(annotated, target, selection=PRIMARY_LAYER)

    restored = read_parsed_levels(target)
    assert list(result.layers) == ["Intensity"]
    assert set(restored.levels["ion"].varm) == {"proteobench:Intensity"}
    record = _object(restored.levels["ion"].metadata["proteobench"])
    provenance = _object(
        _object(_object(restored.metadata["proteobench"])["provenance"])["scoring"]
    )
    assert provenance["selection_mode"] == "primary"
    assert list(_object(record["scoring"])) == ["Intensity"]


def test_named_selection_scores_one_abundance_layer(tmp_path: Path) -> None:
    annotated = _write_annotated_result(tmp_path, _multi_layer_result())
    target = tmp_path / "scored.parquet"

    result = score_result(
        annotated,
        target,
        selection=NamedAbundanceLayer("LFQ/Intensity"),
    )

    restored = read_parsed_levels(target)
    assert list(result.layers) == ["LFQ/Intensity"]
    assert set(restored.levels["ion"].varm) == {"proteobench:LFQ/Intensity"}
    record = _object(restored.levels["ion"].metadata["proteobench"])
    provenance = _object(
        _object(_object(restored.metadata["proteobench"])["provenance"])["scoring"]
    )
    assert provenance["selection_mode"] == "named_abundance"
    assert list(_object(record["scoring"])) == ["LFQ%2FIntensity"]


@pytest.mark.parametrize("suffix", [".h5ad", ".h5mu", ".parquet", ".duckdb"])
def test_all_abundance_layers_round_trip_in_declared_order(
    suffix: str,
    tmp_path: Path,
) -> None:
    annotated = _write_annotated_result(tmp_path, _multi_layer_result(), suffix=suffix)
    target = tmp_path / f"scored{suffix}"

    result = score_result(annotated, target, selection=ALL_ABUNDANCE_LAYERS)

    restored = read_parsed_levels(target)
    assert list(result.layers) == ["Intensity", "LFQ/Intensity"]
    assert list(restored.levels["ion"].varm) == [
        "proteobench:Intensity",
        "proteobench:LFQ/Intensity",
    ]
    record = _object(restored.levels["ion"].metadata["proteobench"])
    layer_records = _object(record["scoring"])
    assert list(layer_records) == ["Intensity", "LFQ%2FIntensity"]
    assert [_object(value)["layer_name"] for value in layer_records.values()] == [
        "Intensity",
        "LFQ/Intensity",
    ]
    assert "QValue" not in layer_records


@pytest.mark.parametrize("name", ["QValue", "Missing"])
def test_named_selection_rejects_layers_without_abundance_role(
    name: str,
    tmp_path: Path,
) -> None:
    annotated = _write_annotated_result(tmp_path, _multi_layer_result())

    with pytest.raises(ValueError, match=r"layer|abundance"):
        score_result(
            annotated,
            tmp_path / "scored.parquet",
            selection=NamedAbundanceLayer(name),
        )


def test_all_abundance_layers_fall_back_to_primary_without_roles(tmp_path: Path) -> None:
    parsed = _multi_layer_result()
    parsed.levels["ion"].uns.pop("layer_roles")
    annotated = _write_annotated_result(tmp_path, parsed)
    target = tmp_path / "scored.parquet"

    result = score_result(annotated, target, selection=ALL_ABUNDANCE_LAYERS)

    assert list(result.layers) == ["Intensity"]
    assert result.selection.fallback == "primary_missing_abundance_roles"
    restored = read_parsed_levels(target)
    provenance = _object(
        _object(_object(restored.metadata["proteobench"])["provenance"])["scoring"]
    )
    assert provenance["selection_fallback"] == "primary_missing_abundance_roles"


@pytest.mark.parametrize(
    ("abundance", "message"),
    [
        (["Intensity", "Intensity"], "duplicate"),
        (["Intensity", "Missing"], "references missing"),
        ("Intensity", "string list"),
    ],
)
def test_all_abundance_layers_reject_corrupt_role_metadata(
    abundance: JsonValue,
    message: str,
    tmp_path: Path,
) -> None:
    parsed = _multi_layer_result()
    parsed.levels["ion"].uns["layer_roles"] = {"abundance": abundance}
    annotated = _write_annotated_result(tmp_path, parsed)

    with pytest.raises(ValueError, match=message):
        score_result(
            annotated,
            tmp_path / "scored.parquet",
            selection=ALL_ABUNDANCE_LAYERS,
        )


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


class _FailsAfterFirstDiagnostics:
    def __init__(self, result: IntermediateResult) -> None:
        self.result = result
        self.calls = 0

    def diagnose(
        self,
        level: QuantitativeLevelInput,
        configuration: ModuleSettings,
        /,
    ) -> IntermediateResult:
        del level, configuration
        self.calls += 1
        if self.calls == 2:
            raise RuntimeError("second layer failed")
        return self.result

    def identity(self) -> dict[str, str]:
        return {"name": "failing-diagnostics", "version": "1"}


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


def test_multi_layer_failure_writes_no_partial_target(tmp_path: Path) -> None:
    annotated = _write_annotated_result(tmp_path, _multi_layer_result())
    target = tmp_path / "scored.parquet"
    diagnostics = MixedSpeciesDiagnostics().diagnose(quantitative_input(), module_settings())

    with pytest.raises(RuntimeError, match="second layer failed"):
        score_result(
            annotated,
            target,
            selection=ALL_ABUNDANCE_LAYERS,
            diagnostic_method=_FailsAfterFirstDiagnostics(diagnostics),
        )

    assert not target.exists()
    restored = read_parsed_levels(annotated)
    assert "scoring" not in _object(restored.levels["ion"].metadata["proteobench"])
    assert not any(name.startswith("proteobench:") for name in restored.levels["ion"].varm)


@pytest.mark.parametrize("owner", ["root", "level"])
def test_score_metadata_alone_prevents_overwrite(tmp_path: Path, owner: str) -> None:
    annotated = _write_annotated_result(tmp_path, parsed_result())
    parsed = read_parsed_levels(annotated)
    root = _object(_object(parsed.metadata["proteobench"])["provenance"])
    local = _object(parsed.levels["ion"].metadata["proteobench"])
    section = root if owner == "root" else local
    section["scoring"] = {}
    source = tmp_path / "existing-scores.parquet"
    target = tmp_path / "refused.parquet"
    write_parsed_levels(parsed, source)
    before = deepcopy((parsed.metadata, parsed.levels["ion"].metadata))
    with pytest.raises(ValueError, match="refusing to overwrite"):
        score_result(source, target)
    assert not target.exists()
    reread = read_parsed_levels(source)
    assert (reread.metadata, reread.levels["ion"].metadata) == before


def test_cli_benchmark_exports_and_reports_verbose_summary(tmp_path: Path) -> None:
    source = tmp_path / "source.parquet"
    benchmarked = tmp_path / "benchmarked.parquet"
    benchmark_export = tmp_path / "benchmark-report" / "result_performance.csv"
    module = tmp_path / "module.toml"
    write_module(module)
    write_parsed_levels(parsed_result(), source)
    messages: list[str] = []
    sink = logger.add(messages.append, format="{message}")
    try:
        assert (
            app(
                [
                    "benchmark",
                    str(source),
                    str(module),
                    str(benchmarked),
                    "--x",
                    "--result-performance",
                    str(benchmark_export),
                    "--verbose",
                ],
                exit_on_error=False,
                result_action="return_value",
            )
            == 0
        )
    finally:
        logger.remove(sink)

    rendered = "".join(messages)
    assert "precursor ion" in benchmark_export.read_text(encoding="utf-8").splitlines()[0]
    proteobot_paths = list(benchmark_export.parent.glob("*.json"))
    assert len(proteobot_paths) == 1
    proteobot = json.loads(proteobot_paths[0].read_text(encoding="utf-8"))
    assert proteobot["software_name"] == "unknown"
    assert proteobot_paths[0].stem == proteobot["intermediate_hash"]
    assert set(proteobot["results"]) == {"1", "2", "3", "4", "5", "6"}
    assert "wrote pMultiQC input" in rendered
    assert "wrote ProteoBot datapoint" in rendered
    assert "ProteoBench score summary" in rendered
    assert "species=['YEAST', 'ECOLI', 'HUMAN']" in rendered
    assert "features total=6 included=3 excluded=3" in rendered
    assert "varm['proteobench:Intensity']" in rendered
    assert "level=ion layer=Intensity" in rendered


def test_cli_default_scores_every_abundance_layer(tmp_path: Path) -> None:
    source = tmp_path / "source.parquet"
    module = tmp_path / "module.toml"
    target = tmp_path / "scored.parquet"
    verbose_target = tmp_path / "scored-verbose.parquet"
    write_module(module)
    write_parsed_levels(_multi_layer_result(), source)
    concise_messages: list[str] = []
    concise_sink = logger.add(concise_messages.append, format="{message}")
    try:
        status = app(
            ["benchmark", str(source), str(module), str(target)],
            exit_on_error=False,
            result_action="return_value",
        )
    finally:
        logger.remove(concise_sink)

    verbose_messages: list[str] = []
    verbose_sink = logger.add(verbose_messages.append, format="{message}")
    try:
        verbose_status = app(
            [
                "benchmark",
                str(source),
                str(module),
                str(verbose_target),
                "--verbose",
            ],
            exit_on_error=False,
            result_action="return_value",
        )
    finally:
        logger.remove(verbose_sink)

    assert status == 0
    assert verbose_status == 0
    assert list(read_parsed_levels(target).levels["ion"].varm) == [
        "proteobench:Intensity",
        "proteobench:LFQ/Intensity",
    ]
    concise = "".join(concise_messages)
    verbose = "".join(verbose_messages)
    for layer_name in ("Intensity", "LFQ/Intensity"):
        assert f"scored level=ion layer={layer_name}" in concise
        assert f"level=ion layer={layer_name}" in verbose


def test_cli_x_scores_only_primary_layer(tmp_path: Path) -> None:
    source = tmp_path / "source.parquet"
    module = tmp_path / "module.toml"
    target = tmp_path / "scored.parquet"
    write_module(module)
    write_parsed_levels(_multi_layer_result(), source)

    status = app(
        ["benchmark", str(source), str(module), str(target), "--x"],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    assert list(read_parsed_levels(target).levels["ion"].varm) == ["proteobench:Intensity"]


@pytest.mark.parametrize(
    "arguments",
    [
        ["benchmark", "source.parquet", "module.toml", "target.parquet"],
        [
            "run",
            "vendor.tsv",
            "proteins.fasta",
            "--params",
            "params.txt",
            "--module",
            "module.toml",
            "--output",
            "target.parquet",
        ],
    ],
)
def test_cli_rejects_named_layer_combined_with_x(arguments: list[str]) -> None:
    status = app(
        [*arguments, "--layer", "Intensity", "--x"],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 2


@pytest.mark.parametrize(
    "arguments",
    [
        ["benchmark", "source.parquet", "module.toml", "target.parquet"],
        [
            "run",
            "vendor.tsv",
            "proteins.fasta",
            "--params",
            "params.txt",
            "--module",
            "module.toml",
            "--output",
            "target.parquet",
        ],
    ],
)
def test_cli_requires_one_layer_for_result_performance(arguments: list[str]) -> None:
    status = app(
        [
            *arguments,
            "--result-performance",
            "result_performance.csv",
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 2


def test_cli_rejects_result_performance_for_non_ion_level(tmp_path: Path) -> None:
    parsed = parsed_result()
    level = parsed.levels.pop("ion")
    level.uns["quantification_level"] = "peptidoform"
    parsed.levels["peptidoform"] = level
    source = tmp_path / "source.parquet"
    module = tmp_path / "module.toml"
    target = tmp_path / "scored.parquet"
    export = tmp_path / "result_performance.csv"
    write_module(module)
    module.write_text(
        module.read_text(encoding="utf-8").replace('level = "ion"', 'level = "peptidoform"'),
        encoding="utf-8",
    )
    write_parsed_levels(parsed, source)

    status = app(
        [
            "benchmark",
            str(source),
            str(module),
            str(target),
            "--x",
            "--result-performance",
            str(export),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 1
    assert not export.exists()


def test_cli_exposes_only_complete_workflows() -> None:
    command_names = {name for name in app.resolved_commands() if not name.startswith("-")}

    assert command_names == {"benchmark", "run"}


@pytest.mark.parametrize(
    ("command", "descriptions"),
    [
        (
            "benchmark",
            (
                "Existing APB2 result to annotate and score",
                "ProteoBench module settings TOML",
                "New scored APB2 result; format selected from its suffix",
                "One named APB abundance layer to score",
                "Score only the APB primary layer represented by X",
                "Write result_performance.csv and sibling",
                "ProteoBot",
                "Report detailed diagnostics and scores",
            ),
        ),
        (
            "run",
            (
                "Vendor result table or directory",
                "One or more protein FASTA files",
                "Vendor search-parameter file",
                "ProteoBench module settings TOML",
                "New scored APB2 .h5ad, .h5mu, .parquet, or",
                ".duckdb result",
                "Vendor software selecting APB2 parsing rules",
                "Software parser override for the parameter file",
                "FASTA peptide-matching backend",
                "Treat isoleucine and leucine as equivalent",
                "Separator between protein accessions",
                "One quantification level to convert",
                "One named APB abundance layer to score",
                "Score only the APB primary layer represented by X",
                "Write result_performance.csv and sibling",
                "ProteoBot JSON",
                "Promote APB2 layer-contract warnings to errors",
                "Report detailed conversion and scoring",
            ),
        ),
    ],
)
def test_cli_help_describes_every_argument_and_parameter(
    command: str,
    descriptions: tuple[str, ...],
    capsys: pytest.CaptureFixture[str],
) -> None:
    app(
        [command, "--help"],
        exit_on_error=False,
        result_action="return_value",
    )

    rendered = capsys.readouterr().out
    assert all(description in rendered for description in descriptions)
    assert "--no-" not in rendered


def _object(value: JsonValue) -> dict[str, JsonValue]:
    assert isinstance(value, dict)
    return value

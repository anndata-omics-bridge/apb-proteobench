"""ProteoBench annotation, APB persistence, protocol, and CLI behavior."""

from __future__ import annotations

import ast
import json
from pathlib import Path

import polars as pl
import pytest
from apb2.api import (
    FinalLayerTable,
    JsonValue,
    ParsedLevels,
    read_parsed_levels,
    write_parsed_levels,
)
from loguru import logger

from apb_proteobench import api as public_api
from apb_proteobench.annotation import ProteoBenchAnnotationParser
from apb_proteobench.api import ProteoBenchAnalysisResult, ProteoBenchAnalyzer
from apb_proteobench.calculation.contracts import QuantitativeLevelInput
from apb_proteobench.calculation.intermediate import IntermediateResult
from apb_proteobench.calculation.metrics import ProteoBenchScores
from apb_proteobench.cli.app import app
from apb_proteobench.configuration.load import load_module
from apb_proteobench.configuration.schema import ModuleSettings
from apb_proteobench.workflow import (
    MixedSpeciesDiagnostics,
    ProteoBenchCompatibleScoring,
    analyze_level,
)
from conftest import module_settings, parsed_result, quantitative_input, write_module


def test_public_api_is_an_in_memory_proteobench_boundary() -> None:
    path = Path(__file__).parents[1] / "src/apb_proteobench/api.py"
    document = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imported = {
        (node.module, name.name)
        for node in ast.walk(document)
        if isinstance(node, ast.ImportFrom) and node.module is not None
        for name in node.names
    }

    assert {item for item in imported if item[0].startswith("apb2")} == {
        ("apb2.api", "ParsedLevels"),
    }
    assert not any(
        module == "pathlib" or module.startswith(("apb_fasta", "protein_fasta"))
        for module, _name in imported
    )


def test_annotation_requires_exact_coverage_and_embeds_complete_configuration(
    tmp_path: Path,
) -> None:
    module = tmp_path / "module.toml"
    write_module(module, alias=True)
    parsed = parsed_result()

    result = ProteoBenchAnnotationParser.from_path(module).parse(parsed).annotate()

    level = result.levels["ion"]
    assert level.obs.frame.columns == ["Run", "sample_name", "condition"]
    assert level.obs.frame.get_column("sample_name").to_list() == ["A1", "A2", "B1", "B2"]
    record = _object(result.metadata["proteobench"])
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

    level = result.levels["ion"]
    assert level.obs.frame.get_column("condition").to_list() == ["A", "A", "B", "B"]
    assert level.obs.frame.get_column("sample_name").to_list() == ["A1", "A2", "B1", "B2"]


def test_annotation_preserves_root_annotation_tables_and_relations(tmp_path: Path) -> None:
    module = tmp_path / "module.toml"
    write_module(module, alias=True)
    parsed = (
        parsed_result()
        .with_annotation_table("members", pl.DataFrame({"member": ["P1"]}), ("member",))
        .with_feature_relation(
            "membership",
            "members",
            "ion",
            pl.DataFrame({"row": [0], "column": [0], "value": [1.0]}),
        )
    )

    result = ProteoBenchAnnotationParser.from_path(module).parse(parsed).annotate()

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

    sdrf = module.with_suffix(".sdrf.tsv")
    missing_sdrf = tmp_path / "missing.sdrf.tsv"
    missing_sdrf.write_text(
        "".join(sdrf.read_text(encoding="utf-8").splitlines(keepends=True)[:-1]),
        encoding="utf-8",
    )
    missing_module = tmp_path / "missing.toml"
    missing_module.write_text(
        module.read_text(encoding="utf-8").replace(sdrf.name, missing_sdrf.name),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="complete sample annotation required"):
        ProteoBenchAnnotationParser.from_path(missing_module).parse(parsed_result())


@pytest.mark.parametrize("suffix", [".h5ad", ".h5mu", ".parquet", ".duckdb"])
def test_annotation_scoring_and_roundtrip_through_every_apb_format(
    suffix: str,
    tmp_path: Path,
) -> None:
    scored = tmp_path / f"scored{suffix}"
    module = tmp_path / "module.toml"
    write_module(module)
    parsed = parsed_result()
    loaded = load_module(module)
    baseline = ProteoBenchAnnotationParser(loaded).parse(parsed).annotate()
    result = ProteoBenchAnalyzer(loaded).analyze(parsed)
    write_parsed_levels(result.parsed, scored)
    restored = read_parsed_levels(scored)

    assert result.layers["Intensity"].analysis.scores.nr_feature == 3
    assert restored.levels["ion"].varm["proteobench:Intensity"].height == 6
    stored = restored.levels["ion"].metadata["proteobench"]
    assert isinstance(stored, dict)
    provenance = _object(_object(restored.metadata["proteobench"])["provenance"])
    assert _object(provenance["scoring"])["schema_version"] == "3"
    assert set(provenance) == {"annotation", "scoring"}
    assert (
        _object(restored.levels["ion"].metadata["prolfquapp"])["annotation"]
        == _object(baseline.levels["ion"].metadata["prolfquapp"])["annotation"]
    )
    assert set(stored) == {"scoring"}
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
    with pytest.raises(ValueError, match="annotation columns already present"):
        ProteoBenchAnalyzer(loaded).analyze(result.parsed)


def test_analyzer_annotates_and_scores_without_physical_io(tmp_path: Path) -> None:
    target = tmp_path / "benchmarked.parquet"
    module = tmp_path / "module.toml"
    write_module(module)

    result = ProteoBenchAnalyzer(load_module(module)).analyze(parsed_result())

    assert result.layers["Intensity"].analysis.scores.nr_feature == 3
    assert "sample_name" in result.parsed.levels["ion"].obs.frame
    assert "proteobench:Intensity" in result.parsed.levels["ion"].varm
    assert not target.exists()


def _multi_layer_result() -> ParsedLevels:
    parsed = parsed_result()
    level = parsed.levels["ion"]
    values = level.layers["Intensity"].values
    level.layers["LFQ/Intensity"] = FinalLayerTable(
        layer_name="LFQ/Intensity",
        values=(values.clone()).drop(("feature",), strict=False),
        semantic_roles=("abundance",),
    )
    level.layers["QValue"] = FinalLayerTable(
        layer_name="QValue",
        values=(values.clone()).drop(("feature",), strict=False),
        semantic_roles=(),
    )
    return parsed


def _analyze_result(
    folder: Path,
    parsed: ParsedLevels,
    /,
    *,
    layers: tuple[str, ...] | None = None,
) -> ProteoBenchAnalysisResult:
    module = folder / "module.toml"
    write_module(module)
    return ProteoBenchAnalyzer(load_module(module), layers=layers).analyze(parsed)


def test_default_selection_scores_every_abundance_layer(tmp_path: Path) -> None:
    target = tmp_path / "scored.parquet"

    result = _analyze_result(tmp_path, _multi_layer_result())
    write_parsed_levels(result.parsed, target)

    restored = read_parsed_levels(target)
    assert list(result.layers) == ["Intensity", "LFQ/Intensity"]
    assert set(restored.levels["ion"].varm) == {
        "fasta_validation",
        "proteobench:Intensity",
        "proteobench:LFQ/Intensity",
    }
    record = _object(restored.levels["ion"].metadata["proteobench"])
    provenance = _object(
        _object(_object(restored.metadata["proteobench"])["provenance"])["scoring"]
    )
    assert provenance["layers"] == ["Intensity", "LFQ/Intensity"]
    assert "resolved_layers" not in provenance
    assert list(_object(record["scoring"])) == ["Intensity", "LFQ%2FIntensity"]


def test_primary_selection_scores_only_x_layer(tmp_path: Path) -> None:
    target = tmp_path / "scored.parquet"

    result = _analyze_result(tmp_path, _multi_layer_result(), layers=("Intensity",))
    write_parsed_levels(result.parsed, target)

    restored = read_parsed_levels(target)
    assert list(result.layers) == ["Intensity"]
    assert set(restored.levels["ion"].varm) == {"fasta_validation", "proteobench:Intensity"}
    record = _object(restored.levels["ion"].metadata["proteobench"])
    provenance = _object(
        _object(_object(restored.metadata["proteobench"])["provenance"])["scoring"]
    )
    assert provenance["layers"] == ["Intensity"]
    assert list(_object(record["scoring"])) == ["Intensity"]


def test_named_selection_scores_one_abundance_layer(tmp_path: Path) -> None:
    target = tmp_path / "scored.parquet"

    result = _analyze_result(
        tmp_path,
        _multi_layer_result(),
        layers=("LFQ/Intensity",),
    )
    write_parsed_levels(result.parsed, target)

    restored = read_parsed_levels(target)
    assert list(result.layers) == ["LFQ/Intensity"]
    assert set(restored.levels["ion"].varm) == {"fasta_validation", "proteobench:LFQ/Intensity"}
    record = _object(restored.levels["ion"].metadata["proteobench"])
    provenance = _object(
        _object(_object(restored.metadata["proteobench"])["provenance"])["scoring"]
    )
    assert provenance["layers"] == ["LFQ/Intensity"]
    assert list(_object(record["scoring"])) == ["LFQ%2FIntensity"]


@pytest.mark.parametrize("suffix", [".h5ad", ".h5mu", ".parquet", ".duckdb"])
def test_all_abundance_layers_round_trip_in_declared_order(
    suffix: str,
    tmp_path: Path,
) -> None:
    target = tmp_path / f"scored{suffix}"

    result = _analyze_result(
        tmp_path,
        _multi_layer_result(),
        layers=None,
    )
    write_parsed_levels(result.parsed, target)

    restored = read_parsed_levels(target)
    assert list(result.layers) == ["Intensity", "LFQ/Intensity"]
    assert list(restored.levels["ion"].varm) == [
        "fasta_validation",
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
    with pytest.raises(ValueError, match=r"layer|abundance"):
        _analyze_result(
            tmp_path,
            _multi_layer_result(),
            layers=(name,),
        )


def test_all_abundance_layers_reject_missing_roles(tmp_path: Path) -> None:
    parsed = _multi_layer_result()
    for layer in parsed.levels["ion"].layers.values():
        layer.semantic_roles = ()
    with pytest.raises(ValueError, match="abundance role"):
        _analyze_result(tmp_path, parsed, layers=None)


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


def test_analyzer_binds_supplied_calculation_methods(tmp_path: Path) -> None:
    module = tmp_path / "module.toml"
    write_module(module)
    configuration = module_settings()
    expected = analyze_level(
        quantitative_input(),
        configuration,
        MixedSpeciesDiagnostics(),
        ProteoBenchCompatibleScoring(),
    )
    diagnostics = _ObservedDiagnostics(expected.diagnostics)
    scoring = _ObservedScoring(expected.scores)
    analyzer = ProteoBenchAnalyzer(
        load_module(module),
        diagnostic_method=diagnostics,
        scoring_method=scoring,
    )

    result = analyzer.analyze(parsed_result())

    assert diagnostics.called and scoring.called
    assert result.layers["Intensity"].analysis.diagnostic_method["name"] == "test-diagnostics"
    assert result.layers["Intensity"].analysis.scoring_method["name"] == "test-scoring"


def test_multi_layer_failure_mutates_no_input_or_target(tmp_path: Path) -> None:
    parsed = _multi_layer_result()
    target = tmp_path / "scored.parquet"
    module = tmp_path / "module.toml"
    write_module(module)
    diagnostics = MixedSpeciesDiagnostics().diagnose(quantitative_input(), module_settings())
    analyzer = ProteoBenchAnalyzer(
        load_module(module),
        layers=None,
        diagnostic_method=_FailsAfterFirstDiagnostics(diagnostics),
    )

    with pytest.raises(RuntimeError, match="second layer failed"):
        analyzer.analyze(parsed)

    assert not target.exists()
    assert "proteobench" not in parsed.metadata
    assert "proteobench" not in parsed.levels["ion"].metadata
    assert not any(name.startswith("proteobench:") for name in parsed.levels["ion"].varm)


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


@pytest.mark.parametrize(
    ("layer_arguments", "expected"),
    [
        ([], ["proteobench:Intensity"]),
        (["--layer", "X"], ["proteobench:Intensity"]),
        (["--layer", "LFQ/Intensity"], ["proteobench:LFQ/Intensity"]),
    ],
)
def test_cli_scores_one_layer_and_defaults_to_x(
    tmp_path: Path,
    layer_arguments: list[str],
    expected: list[str],
) -> None:
    source = tmp_path / "source.parquet"
    module = tmp_path / "module.toml"
    target = tmp_path / "scored.parquet"
    write_module(module)
    write_parsed_levels(_multi_layer_result(), source)
    messages: list[str] = []
    sink = logger.add(messages.append, format="{message}")
    try:
        status = app(
            ["benchmark", str(source), str(module), str(target), *layer_arguments],
            exit_on_error=False,
            result_action="return_value",
        )
    finally:
        logger.remove(sink)

    assert status == 0
    assert list(read_parsed_levels(target).levels["ion"].varm) == ["fasta_validation", *expected]
    assert f"scored level=ion layer={expected[0].removeprefix('proteobench:')}" in "".join(messages)


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
            "--result-performance",
            str(export),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 1
    assert not export.exists()


def test_cli_benchmark_resolves_packaged_module_names(tmp_path: Path) -> None:
    source = tmp_path / "source.parquet"
    write_parsed_levels(parsed_result(), source)
    messages: list[str] = []
    sink = logger.add(messages.append, format="{message}")
    try:
        statuses = [
            app(
                ["benchmark", str(source), name, str(tmp_path / f"{name}.parquet")],
                exit_on_error=False,
                result_action="return_value",
            )
            for name in ("missing", "dda_qexactive")
        ]
    finally:
        logger.remove(sink)

    assert statuses == [1, 1]
    assert "unknown packaged ProteoBench module 'missing'" in messages[0]
    assert "sample annotation matched no observations" in messages[1]


def test_cli_exposes_only_complete_workflows() -> None:
    command_names = {name for name in app.resolved_commands() if not name.startswith("-")}

    assert command_names == {"benchmark", "run"}


@pytest.mark.parametrize(
    ("command", "descriptions"),
    [
        (
            ("benchmark",),
            (
                "Existing APB2 result to annotate and score",
                "Packaged ProteoBench module",
                "New scored APB2 result; format selected from its suffix",
                "Abundance layer to score",
                "Write result_performance.csv and sibling",
                "ProteoBot",
                "Report detailed diagnostics and scores",
            ),
        ),
        (
            ("run", "quant"),
            (
                "Vendor result table or directory",
                "One or more protein FASTA files",
                "Vendor search-parameter file",
                "Packaged ProteoBench module",
                "New scored APB2 .h5ad, .h5mu, .parquet, or",
                ".duckdb result",
                "Parameter-file software; restrict result",
                "recognition to plausible vendors",
                "FASTA peptide-matching backend",
                "Treat isoleucine and leucine as equivalent",
                "Separator between protein accessions",
                "One quantification level to convert",
                "Abundance layer to score",
                "Write result_performance.csv and sibling",
                "ProteoBot JSON",
                "Promote APB2 layer-contract warnings to errors",
                "Report detailed conversion and scoring",
            ),
        ),
        (
            ("run", "entrapment"),
            (
                "Vendor result table or directory",
                "One or more protein FASTA files",
                "Vendor search-parameter file",
                "Packaged entrapment module name",
                "FASTA peptide-matching backend",
                "Promote APB2 layer-contract warnings to errors",
            ),
        ),
        (
            ("run",),
            (
                "quant",
                "Score one quantitative module",
                "entrapment",
                "Score one entrapment module",
            ),
        ),
    ],
)
def test_cli_help_describes_every_argument_and_parameter(
    command: tuple[str, ...],
    descriptions: tuple[str, ...],
    capsys: pytest.CaptureFixture[str],
) -> None:
    app(
        [*command, "--help"],
        exit_on_error=False,
        result_action="return_value",
    )

    rendered = capsys.readouterr().out
    assert all(description in rendered for description in descriptions)
    assert "--no-" not in rendered
    assert "--params-software" not in rendered
    assert "--x " not in rendered


def _object(value: JsonValue) -> dict[str, JsonValue]:
    assert isinstance(value, dict)
    return value


def test_api_exports_exactly_the_approved_names() -> None:
    assert sorted(public_api.__all__) == [
        "EntrapmentAnalysisResult",
        "EntrapmentAnalyzer",
        "LoadedModule",
        "ProteoBenchAnalysisResult",
        "ProteoBenchAnalyzer",
        "SubmissionContent",
        "load_module",
        "load_packaged_entrapment_module",
        "load_packaged_module",
    ]

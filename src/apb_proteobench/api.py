"""Vendor conversion and result-I/O workflows exposed to Python and the CLI."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

from apb2.annotation_extension import AnnotationResult
from apb2.api import (
    ParsedLevels,
    ParseRuleCompiler,
    QuantificationLevel,
    read_parsed_levels,
    write_parsed_levels,
)
from apb2.result_facade import JsonValue
from apb_fasta.annotation import FastaAnnotationParser
from apb_fasta.calculation.results import FastaAnnotationReports
from apb_fasta.configuration import (
    DEFAULT_FASTA_ANNOTATION_PARAMETERS,
    FastaAnnotationParameters,
)
from protein_fasta.frame import ProteinDatabase, refseq, uniprotkb

from apb_proteobench.annotation import ProteoBenchAnnotationParser
from apb_proteobench.configuration.schema import ModuleSettings
from apb_proteobench.integration import (
    ALL_ABUNDANCE_LAYERS,
    LayerSelection,
    ResolvedLayerSelection,
    ScoredLayerResult,
    diagnostics_slot,
    embedded_configuration,
    extract_layer,
    persist_results,
    resolve_layer_selection,
)
from apb_proteobench.workflow import (
    DiagnosticMethod,
    MixedSpeciesDiagnostics,
    ProteoBenchCompatibleScoring,
    ScoringMethod,
    analyze_level,
)

_DEFAULT_DIAGNOSTICS = MixedSpeciesDiagnostics()
_DEFAULT_SCORING = ProteoBenchCompatibleScoring()
type ValidationChecks = Literal["standard", "strict"]


@dataclass(frozen=True, slots=True)
class ConvertedVendorResult:
    """A vendor conversion and the in-memory APB2 result it produced."""

    input_path: Path
    parameters_path: Path
    output_path: Path
    software: str
    software_version: str | None
    parsed: ParsedLevels


@dataclass(frozen=True, slots=True)
class ScoredResult:
    """Typed workflow evidence used by presenters and API clients."""

    input_path: Path
    output_path: Path
    configuration: ModuleSettings
    selection: ResolvedLayerSelection
    layers: dict[str, ScoredLayerResult]
    search_parameters: dict[str, JsonValue]


@dataclass(frozen=True, slots=True)
class VendorBenchmarkResult:
    """Complete raw-vendor workflow evidence and the final APB2 result."""

    software: str
    software_version: str | None
    parsed: ParsedLevels
    fasta_reports: FastaAnnotationReports
    scored: ScoredResult


@dataclass(frozen=True, slots=True)
class _AnalyzedResult:
    parsed: ParsedLevels
    configuration: ModuleSettings
    selection: ResolvedLayerSelection
    layers: dict[str, ScoredLayerResult]


def _persist_scored(
    analyzed: _AnalyzedResult,
    source: Path,
    target: Path,
    /,
    *,
    search_parameters: dict[str, JsonValue],
) -> ScoredResult:
    """Write the analyzed result and wrap it as the public scoring evidence."""
    target.parent.mkdir(parents=True, exist_ok=True)
    write_parsed_levels(analyzed.parsed, target)
    return ScoredResult(
        input_path=source,
        output_path=target,
        configuration=analyzed.configuration,
        selection=analyzed.selection,
        layers=analyzed.layers,
        search_parameters=search_parameters,
    )


def convert_vendor_result(
    data: Path,
    parameters_path: Path,
    target: Path,
    /,
    *,
    level: QuantificationLevel | None = None,
    software: str | None = None,
    parameters_software: str | None = None,
    checks: ValidationChecks = "standard",
) -> ConvertedVendorResult:
    """Parse a vendor table with APB2's compiler and persist one APB2 result.

    Args:
        data: Vendor result table or canonical multi-file result directory.
        parameters_path: Vendor search-parameter file.
        target: Exact output path in a format supporting the requested level count.
        level: One quantification level. If omitted, parse every compatible level.
        software: Optional vendor slug used to select and verify the packaged rule document.
        parameters_software: Optional independent parameter-parser slug.
        checks: AnnData layer-contract validation level.

    Returns:
        The detected software metadata and parsed in-memory APB2 result.

    Raises:
        ValueError: The inputs cannot select or satisfy one packaged rule document, or the output
            format cannot store the requested conversion.
    """
    _require_new_target(data, target)
    compiler = ParseRuleCompiler(
        data,
        parameters_path,
        requested_levels=None if level is None else (level,),
        checks=checks,
        software=software,
        parameters_software=parameters_software,
    )
    parsed = compiler.compile().parse()
    target.parent.mkdir(parents=True, exist_ok=True)
    write_parsed_levels(parsed, target)
    return ConvertedVendorResult(
        input_path=data,
        parameters_path=parameters_path,
        output_path=target,
        software=compiler.detection.software,
        software_version=compiler.detection.version,
        parsed=parsed,
    )


def annotate_result(source: Path, module: Path, target: Path, /) -> AnnotationResult:
    """Bind module samples/configuration to an APB result and persist a new result."""
    _require_new_target(source, target)
    parsed = read_parsed_levels(source)
    annotation = ProteoBenchAnnotationParser.from_path(module).parse(parsed)
    result = annotation.annotate()
    target.parent.mkdir(parents=True, exist_ok=True)
    write_parsed_levels(result.parsed, target)
    return result


def benchmark_result(
    source: Path,
    module: Path,
    target: Path,
    /,
    *,
    selection: LayerSelection = ALL_ABUNDANCE_LAYERS,
    diagnostic_method: DiagnosticMethod = _DEFAULT_DIAGNOSTICS,
    scoring_method: ScoringMethod = _DEFAULT_SCORING,
) -> ScoredResult:
    """Annotate and score selected layers in one existing APB2 result."""
    _require_new_target(source, target)
    parsed = read_parsed_levels(source)
    annotated = ProteoBenchAnnotationParser.from_path(module).parse(parsed).annotate()
    analyzed = _analyze_parsed(
        annotated.parsed,
        selection=selection,
        diagnostic_method=diagnostic_method,
        scoring_method=scoring_method,
    )
    return _persist_scored(analyzed, source, target, search_parameters={})


def score_result(
    source: Path,
    target: Path,
    /,
    *,
    selection: LayerSelection = ALL_ABUNDANCE_LAYERS,
    diagnostic_method: DiagnosticMethod = _DEFAULT_DIAGNOSTICS,
    scoring_method: ScoringMethod = _DEFAULT_SCORING,
) -> ScoredResult:
    """Score selected layers in the configured APB level and persist a new result."""
    _require_new_target(source, target)
    parsed = read_parsed_levels(source)
    analyzed = _analyze_parsed(
        parsed,
        selection=selection,
        diagnostic_method=diagnostic_method,
        scoring_method=scoring_method,
    )
    return _persist_scored(analyzed, source, target, search_parameters={})


def run_vendor_benchmark(
    data: Path,
    parameters_path: Path,
    fasta_paths: tuple[Path, ...],
    module: Path,
    target: Path,
    /,
    *,
    level: QuantificationLevel | None = None,
    software: str | None = None,
    parameters_software: str | None = None,
    checks: ValidationChecks = "standard",
    selection: LayerSelection = ALL_ABUNDANCE_LAYERS,
    fasta_parameters: FastaAnnotationParameters = DEFAULT_FASTA_ANNOTATION_PARAMETERS,
    diagnostic_method: DiagnosticMethod = _DEFAULT_DIAGNOSTICS,
    scoring_method: ScoringMethod = _DEFAULT_SCORING,
) -> VendorBenchmarkResult:
    """Convert, verify, annotate, score, and persist one APB2 result.

    Quantitative aggregation is deliberately not part of this pipeline. Run the separate
    ``apb-aggregate`` command between conversion and benchmarking when a scored level must
    be derived from a lower one.

    Args:
        data: Vendor result table or canonical multi-file result directory.
        parameters_path: Vendor search-parameter file.
        fasta_paths: One or more protein FASTA files.
        module: ProteoBench module settings with the complete sample design.
        target: Exact final output path in a format supporting the requested level count.
        level: One quantification level. If omitted, parse every compatible level.
        software: Optional vendor slug used to verify packaged-rule detection.
        parameters_software: Optional independent parameter-parser slug.
        checks: AnnData layer-contract validation level.
        selection: Layer-selection policy; defaults to every declared abundance layer.
        fasta_parameters: Peptide matching and reported-assignment settings.
        diagnostic_method: ProteoBench diagnostic implementation.
        scoring_method: ProteoBench scoring implementation.

    Returns:
        The final APB2 value plus conversion, FASTA, and scoring evidence.

    Raises:
        ValueError: Inputs are incomplete, incompatible, or target an unsafe output.
    """
    _require_new_target(data, target)
    if not fasta_paths:
        raise ValueError("at least one FASTA path is required")
    compiler = ParseRuleCompiler(
        data,
        parameters_path,
        requested_levels=None if level is None else (level,),
        checks=checks,
        software=software,
        parameters_software=parameters_software,
    )
    parsed = compiler.compile().parse()
    proteins = ProteinDatabase(uniprotkb, refseq).parse(fasta_paths)
    verified = FastaAnnotationParser(
        parsed,
        proteins,
        parameters=fasta_parameters,
    ).verify_peptides()
    fasta_reports = verified.reports
    parsed = verified.parsed
    parsed = ProteoBenchAnnotationParser.from_path(module).parse(parsed).annotate().parsed
    analyzed = _analyze_parsed(
        parsed,
        selection=selection,
        diagnostic_method=diagnostic_method,
        scoring_method=scoring_method,
    )
    return VendorBenchmarkResult(
        software=compiler.detection.software,
        software_version=compiler.detection.version,
        parsed=analyzed.parsed,
        fasta_reports=fasta_reports,
        scored=_persist_scored(
            analyzed,
            data,
            target,
            search_parameters=cast(
                dict[str, JsonValue],
                compiler.parameters.model_dump(mode="json"),
            ),
        ),
    )


def _analyze_parsed(
    parsed: ParsedLevels,
    /,
    *,
    selection: LayerSelection,
    diagnostic_method: DiagnosticMethod,
    scoring_method: ScoringMethod,
) -> _AnalyzedResult:
    configuration = embedded_configuration(parsed)
    resolved = resolve_layer_selection(parsed, configuration, selection)
    layers: dict[str, ScoredLayerResult] = {}
    for layer_name in resolved.layer_names:
        selected = extract_layer(parsed, configuration, layer_name)
        analysis = analyze_level(
            selected.calculation,
            configuration,
            diagnostic_method,
            scoring_method,
        )
        layers[selected.layer_name] = ScoredLayerResult(
            level_name=selected.level_name,
            layer_name=selected.layer_name,
            diagnostics_slot=diagnostics_slot(selected.layer_name),
            roles=selected.roles,
            analysis=analysis,
        )
    return _AnalyzedResult(
        parsed=persist_results(parsed, resolved, layers),
        configuration=configuration,
        selection=resolved,
        layers=layers,
    )


def _require_new_target(source: Path, target: Path) -> None:
    if source.resolve() == target.resolve():
        raise ValueError("output must differ from input")
    if target.exists():
        raise ValueError(f"output already exists: {target}")

"""Minimal Cyclopts interface for complete ProteoBench workflows."""

from __future__ import annotations

import json
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Annotated, Literal

from apb2.api import (
    ParsedLevels,
    ParseRuleCompiler,
    QuantificationLevel,
    read_parsed_levels,
    write_parsed_levels,
)
from apb_fasta.api import (
    FastaAnnotationParameters,
    FastaAnnotationReports,
    FastaAnnotationResult,
    FastaAnnotator,
)
from cyclopts import App, Parameter
from loguru import logger
from pydantic import ValidationError

from apb_proteobench.api import (
    EntrapmentAnalysisResult,
    EntrapmentAnalyzer,
    LoadedModule,
    ProteoBenchAnalysisResult,
    ProteoBenchAnalyzer,
    load_module,
    load_packaged_entrapment_module,
    load_packaged_module,
)
from apb_proteobench.cli.presentation import report_score
from apb_proteobench.cli.result_performance import write_result_performance_bundle
from apb_proteobench.cli.timings import write_tool_timings

PRIMARY_LAYER_NAME = "X"
"""The ``--layer`` value that selects the APB primary layer stored in ``X``."""

app = App(
    name="apb-proteobench",
    help="Run complete ProteoBench workflows with optional pMultiQC/ProteoBot export",
    help_on_error=True,
)


@dataclass(frozen=True, slots=True)
class RunCliOptions:
    """Options for the complete raw-vendor benchmark workflow."""

    params: Annotated[
        Path | None,
        Parameter(help="Vendor search-parameter file"),
    ] = None
    module: Annotated[
        str | None,
        Parameter(help="Packaged ProteoBench module name, or a module TOML path"),
    ] = None
    output: Annotated[
        Path | None,
        Parameter(help="New scored APB2 .h5ad, .h5mu, .parquet, or .duckdb result"),
    ] = None
    software: Annotated[
        str | None,
        Parameter(help="Parameter-file software; restrict result recognition to plausible vendors"),
    ] = None
    backend: Annotated[
        Literal["auto", "ahocorapy", "ahocorasick_rs"],
        Parameter(help="FASTA peptide-matching backend"),
    ] = "auto"
    il_equivalent: Annotated[
        bool,
        Parameter(negative=False, help="Treat isoleucine and leucine as equivalent"),
    ] = False
    protein_group_separator: Annotated[
        str,
        Parameter(help="Separator between protein accessions"),
    ] = ";"
    level: Annotated[
        QuantificationLevel | None,
        Parameter(help="One quantification level to convert"),
    ] = None
    layer: Annotated[
        str,
        Parameter(help='Abundance layer to score; "X" is the APB primary layer'),
    ] = PRIMARY_LAYER_NAME
    result_performance: Annotated[
        Path | None,
        Parameter(help="Write result_performance.csv and sibling ProteoBot JSON"),
    ] = None
    timings_dir: Annotated[
        Path | None,
        Parameter(help="Write separate APB2, FASTA, and ProteoBench timing JSON files"),
    ] = None
    strict: Annotated[
        bool,
        Parameter(negative=False, help="Promote APB2 layer-contract warnings to errors"),
    ] = False


DEFAULT_RUN_CLI_OPTIONS = RunCliOptions()


@dataclass(frozen=True, slots=True)
class EntrapmentCliOptions:
    """Options for the complete raw-vendor entrapment workflow."""

    params: Annotated[
        Path | None,
        Parameter(help="Vendor search-parameter file"),
    ] = None
    module: Annotated[
        str,
        Parameter(help="Packaged entrapment module name"),
    ] = "entrapment_dia_astral"
    output: Annotated[
        Path | None,
        Parameter(help="New scored APB2 .h5ad, .h5mu, .parquet, or .duckdb result"),
    ] = None
    software: Annotated[
        str | None,
        Parameter(help="Parameter-file software; restrict result recognition to plausible vendors"),
    ] = None
    backend: Annotated[
        Literal["auto", "ahocorapy", "ahocorasick_rs"],
        Parameter(help="FASTA peptide-matching backend"),
    ] = "auto"
    il_equivalent: Annotated[
        bool,
        Parameter(negative=False, help="Treat isoleucine and leucine as equivalent"),
    ] = False
    protein_group_separator: Annotated[
        str,
        Parameter(help="Separator between protein accessions"),
    ] = ";"
    timings_dir: Annotated[
        Path | None,
        Parameter(help="Write separate APB2, FASTA, and ProteoBench timing JSON files"),
    ] = None
    strict: Annotated[
        bool,
        Parameter(negative=False, help="Promote APB2 layer-contract warnings to errors"),
    ] = False


DEFAULT_ENTRAPMENT_CLI_OPTIONS = EntrapmentCliOptions()

run_app = App(
    name="run",
    help="Convert vendor files, verify peptides against FASTA, and score one ProteoBench module",
    help_on_error=True,
)
app.command(run_app)


@dataclass(frozen=True, slots=True)
class _VerifiedConversion:
    """Vendor files converted by APB2, their peptides verified against FASTA, and timings."""

    compiler: ParseRuleCompiler
    fasta: FastaAnnotator
    verified: FastaAnnotationResult
    seconds: dict[str, float]
    levels: tuple[dict[str, str | float], ...]


def _convert_and_verify(
    data: Path,
    params: Path,
    fasta_paths: tuple[Path, ...],
    /,
    *,
    level: QuantificationLevel | None,
    software: str | None,
    strict: bool,
    backend: Literal["auto", "ahocorapy", "ahocorasick_rs"],
    il_equivalent: bool,
    protein_group_separator: str,
) -> _VerifiedConversion:
    """Convert vendor files with APB2 and verify their peptides against the FASTA files."""
    if not fasta_paths:
        raise ValueError("at least one FASTA path is required")
    started = perf_counter()
    compiler = ParseRuleCompiler(
        data,
        params,
        requested_levels=None if level is None else (level,),
        checks="strict" if strict else "standard",
        software=software,
    )
    parser = compiler.compile()
    compile_seconds = perf_counter() - started
    started = perf_counter()
    parsed, level_timings = parser.parse_with_timings()
    read_seconds = sum(timing.read_seconds for timing in level_timings)
    parse_seconds = max(0.0, perf_counter() - started - read_seconds)
    started = perf_counter()
    fasta = FastaAnnotator.read(
        fasta_paths,
        parameters=FastaAnnotationParameters(
            protein_group_separator=protein_group_separator,
            matcher_backend=backend,
            il_equivalent=il_equivalent,
        ),
    )
    fasta_load_seconds = perf_counter() - started
    started = perf_counter()
    verified = fasta.verify_peptides(parsed)
    return _VerifiedConversion(
        compiler=compiler,
        fasta=fasta,
        verified=verified,
        seconds={
            "compile": compile_seconds,
            "read": read_seconds,
            "parse": parse_seconds,
            "load_database": fasta_load_seconds,
            "verify_peptides": perf_counter() - started,
        },
        levels=tuple(
            {
                "level": timing.level,
                "read_seconds": timing.read_seconds,
                "parse_seconds": timing.parse_seconds,
            }
            for timing in level_timings
        ),
    )


@app.command
def benchmark(
    source: Annotated[Path, Parameter(help="Existing APB2 result to annotate and score")],
    module: Annotated[
        str,
        Parameter(help="Packaged ProteoBench module name, or a module TOML path"),
    ],
    target: Annotated[
        Path,
        Parameter(help="New scored APB2 result; format selected from its suffix"),
    ],
    /,
    *,
    layer: Annotated[
        str,
        Parameter(help='Abundance layer to score; "X" is the APB primary layer'),
    ] = PRIMARY_LAYER_NAME,
    result_performance: Annotated[
        Path | None,
        Parameter(help="Write result_performance.csv and sibling ProteoBot JSON"),
    ] = None,
    verbose: Annotated[
        bool,
        Parameter(negative=False, help="Report detailed diagnostics and scores"),
    ] = False,
) -> int:
    """Write scored APB2 TARGET and optional pMultiQC/ProteoBot export bundle."""
    try:
        _require_new_target(source, target)
        parsed = read_parsed_levels(source)
        loaded = _load_module(module)
        result = ProteoBenchAnalyzer(loaded, layers=_layer_names(layer, parsed, loaded)).analyze(
            parsed
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        write_parsed_levels(result.parsed, target)
        _export_result_performance(result, {}, result_performance)
    except (OSError, ValueError, ValidationError) as error:
        logger.error(str(error))
        return 1
    report_score(result, source, target, verbose=verbose)
    return 0


@run_app.command(name="quant")
def run_quant(
    data: Annotated[Path, Parameter(help="Vendor result table or directory")],
    /,
    *fasta_paths: Annotated[
        Path,
        Parameter(help="One or more protein FASTA files, or their protein-fasta database Parquet"),
    ],
    options: Annotated[RunCliOptions, Parameter(name="*")] = DEFAULT_RUN_CLI_OPTIONS,
    verbose: Annotated[
        bool,
        Parameter(negative=False, help="Report detailed conversion and scoring diagnostics"),
    ] = False,
) -> int:
    """Score one quantitative module: write scored APB2 output and optional pMultiQC/ProteoBot export.

    Quantitative aggregation is a separate step: run apb-aggregate between conversion
    and benchmarking when the scored level must be derived from a lower one. Use
    --level LEVEL to convert one quantification level. --layer NAME scores one abundance
    layer; the default "X" is the APB primary layer.
    --timings-dir records internal operation timings outside the scored result.
    """
    if options.params is None:
        logger.error("pass --params PATH for the vendor search-parameter file")
        return 1
    if options.module is None:
        logger.error("pass --module NAME or --module PATH for the ProteoBench module")
        return 1
    if options.output is None:
        logger.error("pass --output PATH for the final APB2 result")
        return 1
    try:
        _require_new_target(data, options.output)
        timing_targets = _run_timing_targets(options.timings_dir)
        _require_new_timing_targets(timing_targets)
        conversion = _convert_and_verify(
            data,
            options.params,
            fasta_paths,
            level=options.level,
            software=options.software,
            strict=options.strict,
            backend=options.backend,
            il_equivalent=options.il_equivalent,
            protein_group_separator=options.protein_group_separator,
        )
        seconds = dict(conversion.seconds)
        started = perf_counter()
        loaded_module = _load_module(options.module)
        seconds["load_module"] = perf_counter() - started
        started = perf_counter()
        parsed = conversion.verified.parsed
        result = ProteoBenchAnalyzer(
            loaded_module, layers=_layer_names(options.layer, parsed, loaded_module)
        ).analyze(parsed)
        seconds["analyze"] = perf_counter() - started
        started = perf_counter()
        options.output.parent.mkdir(parents=True, exist_ok=True)
        write_parsed_levels(result.parsed, options.output)
        seconds["write"] = perf_counter() - started
        search_parameters = conversion.compiler.parameters.model_dump(mode="json")
        started = perf_counter()
        _export_result_performance(result, search_parameters, options.result_performance)
        seconds["export"] = perf_counter() - started
        _write_run_timings(
            timing_targets,
            seconds,
            conversion.levels,
            exported=options.result_performance is not None,
        )
    except (OSError, ValueError, ValidationError) as error:
        logger.error(str(error))
        return 1
    _report_vendor_benchmark(
        conversion.compiler.detection.software,
        conversion.compiler.detection.version,
        conversion.verified.reports,
        result,
        data,
        options.output,
        verbose=verbose,
    )
    return 0


@run_app.command(name="entrapment")
def run_entrapment(
    data: Annotated[Path, Parameter(help="Vendor result table or directory")],
    /,
    *fasta_paths: Annotated[
        Path,
        Parameter(help="One or more protein FASTA files, or their protein-fasta database Parquet"),
    ],
    options: Annotated[EntrapmentCliOptions, Parameter(name="*")] = DEFAULT_ENTRAPMENT_CLI_OPTIONS,
) -> int:
    """Score one entrapment module: FDP estimates for every precursor q-value kind.

    The result offers its q-value kinds through apb-catalog's proteobench_entrapment set;
    each kind is scored. Labels and pairs come from the entrapment FASTA, or from its Parquet
    protein database written by protein-fasta database.
    --timings-dir records internal operation timings outside the scored result.
    """
    if options.params is None:
        logger.error("pass --params PATH for the vendor search-parameter file")
        return 1
    if options.output is None:
        logger.error("pass --output PATH for the final APB2 result")
        return 1
    try:
        _require_new_target(data, options.output)
        timing_targets = _run_timing_targets(options.timings_dir)
        _require_new_timing_targets(timing_targets)
        started = perf_counter()
        settings = load_packaged_entrapment_module(options.module)
        module_seconds = perf_counter() - started
        conversion = _convert_and_verify(
            data,
            options.params,
            fasta_paths,
            level=settings.level,
            software=options.software,
            strict=options.strict,
            backend=options.backend,
            il_equivalent=options.il_equivalent,
            protein_group_separator=options.protein_group_separator,
        )
        seconds = dict(conversion.seconds)
        seconds["load_module"] = module_seconds
        started = perf_counter()
        result = EntrapmentAnalyzer(settings, conversion.fasta).analyze(conversion.verified.parsed)
        seconds["analyze"] = perf_counter() - started
        started = perf_counter()
        options.output.parent.mkdir(parents=True, exist_ok=True)
        write_parsed_levels(result.parsed, options.output)
        seconds["write"] = perf_counter() - started
        _write_run_timings(timing_targets, seconds, conversion.levels, exported=False)
    except (OSError, ValueError, ValidationError, LookupError) as error:
        logger.error(str(error))
        return 1
    _report_entrapment(conversion, result, options.output)
    return 0


def _report_entrapment(
    conversion: _VerifiedConversion, result: EntrapmentAnalysisResult, target: Path, /
) -> None:
    logger.info(
        "vendor={} software_version={}",
        conversion.compiler.detection.software,
        conversion.compiler.detection.version or "missing",
    )
    for kind, scores in result.scores.items():
        logger.info(
            "kind={} precursors={} reported_fdr={:.6g} lower_bound_FDP={:.6g} "
            "combined_FDP={:.6g} ({}) paired_FDP={:.6g} ({})",
            kind,
            scores.nr_id_features,
            scores.reported_fdr_parsed_from_input,
            scores.lower_bound_FDP,
            scores.combined_FDP,
            scores.category_combined,
            scores.paired_FDP,
            scores.category_paired,
        )
    logger.info("output={}", target)


def _load_module(module: str, /) -> LoadedModule:
    """Load a module TOML path, or a packaged module when the value is not a TOML path."""
    if module.endswith(".toml"):
        return load_module(Path(module))
    return load_packaged_module(module)


def _run_timing_targets(directory: Path | None) -> dict[str, Path]:
    """Name each optional timing artifact of the integrated run."""
    if directory is None:
        return {}
    return {
        "apb2": directory / "apb2.convert.timings.json",
        "fasta": directory / "apb-fasta.verify-peptides.timings.json",
        "proteobench": directory / "apb-proteobench.benchmark.timings.json",
    }


def _require_new_timing_targets(targets: Mapping[str, Path]) -> None:
    """Refuse a timing collision before creating a scientific result."""
    for target in targets.values():
        if target.exists():
            raise ValueError(f"timing output already exists: {target}")


def _write_run_timings(
    targets: Mapping[str, Path],
    seconds: Mapping[str, float],
    levels: tuple[dict[str, str | float], ...],
    /,
    *,
    exported: bool,
) -> None:
    """Publish one independent timing document for each integrated operation."""
    if not targets:
        return
    write_tool_timings(
        targets["apb2"],
        tool="apb2",
        operation="convert",
        phases=tuple((name, seconds[name]) for name in ("compile", "read", "parse")),
        levels=levels,
    )
    write_tool_timings(
        targets["fasta"],
        tool="apb-fasta",
        operation="verify-peptides",
        phases=tuple((name, seconds[name]) for name in ("load_database", "verify_peptides")),
    )
    names = ["load_module", "analyze", "write"]
    if exported:
        names.append("export")
    write_tool_timings(
        targets["proteobench"],
        tool="apb-proteobench",
        operation="benchmark",
        phases=tuple((name, seconds[name]) for name in names),
    )


def _report_vendor_benchmark(
    software: str,
    software_version: str | None,
    fasta_reports: FastaAnnotationReports,
    result: ProteoBenchAnalysisResult,
    source: Path,
    target: Path,
    /,
    *,
    verbose: bool,
) -> None:
    logger.info(
        "vendor={} software_version={}",
        software,
        software_version or "missing",
    )
    for level, coverage in fasta_reports.peptide_levels.items():
        logger.info(
            "level={} peptides_in_fasta={}/{} unmatched={}",
            level,
            coverage.matched_feature_count,
            coverage.feature_count,
            coverage.unmatched_feature_count,
        )
    report_score(result, source, target, verbose=verbose)


def _layer_names(
    layer: str, parsed: ParsedLevels, module: LoadedModule, /
) -> tuple[str, ...] | None:
    """Name the scored layer; "X" is the configured level's primary layer.

    An absent level returns ``None``: the analyzer then reports the unavailable level.
    """
    if layer != PRIMARY_LAYER_NAME:
        return (layer,)
    level = parsed.levels.get(module.settings.general.level)
    return None if level is None else (level.primary_layer_name,)


def _export_result_performance(
    result: ProteoBenchAnalysisResult,
    search_parameters: Mapping[str, object],
    target: Path | None,
    /,
) -> None:
    if target is None:
        return
    if len(result.layers) != 1:
        raise ValueError("result-performance export requires exactly one selected layer")
    selected = next(iter(result.layers.values()))
    if selected.level_name != "ion":
        raise ValueError(
            "result-performance export currently supports only the ion level; "
            f"got {selected.level_name!r}"
        )
    datapoint = _proteobot_datapoint(result, search_parameters, selected.layer_name)
    written = write_result_performance_bundle(
        selected.analysis.diagnostics.legacy,
        datapoint,
        target,
        content=result.submission(selected.layer_name),
    )
    logger.info("wrote pMultiQC input {}", written.csv)
    logger.info("wrote ProteoBot datapoint {}", written.proteobot_json)


def _proteobot_datapoint(
    result: ProteoBenchAnalysisResult,
    search_parameters: Mapping[str, object],
    layer_name: str,
    /,
) -> dict[str, object]:
    selected = result.layers[layer_name]
    score_document: object = json.loads(selected.analysis.scores.model_dump_json())
    if not isinstance(score_document, dict):
        raise TypeError("ProteoBench scores did not serialize to a JSON object")
    parameters = search_parameters
    software_name = _text_parameter(parameters, "software_name") or "unknown"
    datapoint: dict[str, object] = {
        "software_name": software_name,
        "software_version": _text_parameter(parameters, "software_version") or "",
        "search_engine": _text_parameter(parameters, "search_engine"),
        "search_engine_version": _text_parameter(parameters, "search_engine_version"),
        "ident_fdr_psm": _probability_parameter(parameters, "ident_fdr_psm"),
        "ident_fdr_peptide": _probability_parameter(parameters, "ident_fdr_peptide"),
        "ident_fdr_protein": _probability_parameter(parameters, "ident_fdr_protein"),
        "enable_match_between_runs": _bool_parameter(
            parameters, "enable_match_between_runs", default=False
        ),
        "precursor_mass_tolerance": _tolerance_parameter(parameters, "precursor_mass_tolerance"),
        "fragment_mass_tolerance": _tolerance_parameter(parameters, "fragment_mass_tolerance"),
        "enzyme": _text_parameter(parameters, "enzyme"),
        "allowed_miscleavages": _scalar_parameter(parameters, "allowed_miscleavages"),
        "min_peptide_length": _scalar_parameter(parameters, "min_peptide_length"),
        "max_peptide_length": _scalar_parameter(parameters, "max_peptide_length"),
        "is_temporary": False,
        **score_document,
        "comments": "",
        "old_new": "new",
        "semi_enzymatic": _bool_parameter(parameters, "semi_enzymatic", default=False),
        "fixed_mods": _modifications_parameter(parameters, "fixed_mods"),
        "variable_mods": _modifications_parameter(parameters, "variable_mods"),
        "max_mods": _scalar_parameter(parameters, "max_mods"),
        "min_precursor_charge": _scalar_parameter(parameters, "min_precursor_charge"),
        "max_precursor_charge": _scalar_parameter(parameters, "max_precursor_charge"),
        "quantification_method": _text_parameter(parameters, "quantification_method"),
        "protein_inference": _text_parameter(parameters, "protein_inference"),
        "abundance_normalization_ions": _scalar_parameter(
            parameters, "abundance_normalization_ions"
        ),
        "postprocessing_performed": False,
        "postprocessing_description": None,
    }
    for name in (
        "min_precursor_mz",
        "max_precursor_mz",
        "min_fragment_mz",
        "max_fragment_mz",
        "scan_window",
        "predictors_library",
    ):
        if name in parameters:
            datapoint[name] = _scalar_parameter(parameters, name)
    datapoint["submission_comments"] = ""
    return datapoint


def _scalar_parameter(parameters: Mapping[str, object], name: str, /) -> object:
    return parameters.get(name)


def _text_parameter(parameters: Mapping[str, object], name: str, /) -> str | None:
    value = parameters.get(name)
    return value if isinstance(value, str) and value else None


def _bool_parameter(
    parameters: Mapping[str, object],
    name: str,
    /,
    *,
    default: bool,
) -> bool:
    value = parameters.get(name)
    return value if isinstance(value, bool) else default


def _probability_parameter(parameters: Mapping[str, object], name: str, /) -> object:
    value = parameters.get(name)
    if isinstance(value, dict):
        return value.get("value")
    return value


def _tolerance_parameter(parameters: Mapping[str, object], name: str, /) -> object:
    value = parameters.get(name)
    if not isinstance(value, dict):
        return value
    if value.get("mode") == "automatic":
        return value.get("label") or "Automatic"
    magnitude = value.get("value")
    unit = value.get("unit")
    if not isinstance(magnitude, int | float) or not isinstance(unit, str):
        return None
    rendered = f"{magnitude:g}"
    return f"[-{rendered} {unit}, {rendered} {unit}]"


def _modifications_parameter(parameters: Mapping[str, object], name: str, /) -> str | None:
    value = parameters.get(name)
    if not isinstance(value, list):
        return value if isinstance(value, str) and value else None
    names: list[str] = []
    for modification in value:
        if not isinstance(modification, dict):
            continue
        modification_name = modification.get("name")
        if isinstance(modification_name, str):
            names.append(modification_name)
    return ", ".join(names) if names else None


def _require_new_target(source: Path, target: Path, /) -> None:
    if source.resolve() == target.resolve():
        raise ValueError("output must differ from input")
    if target.exists():
        raise ValueError(f"output already exists: {target}")


def main() -> int:
    """Console-script entry point."""
    result = app()
    return int(result) if result is not None else 0


if __name__ == "__main__":
    sys.exit(main())

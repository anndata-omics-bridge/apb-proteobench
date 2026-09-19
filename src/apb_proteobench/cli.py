"""Minimal Cyclopts interface for complete ProteoBench workflows."""

from __future__ import annotations

import json
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Literal

from apb2.api import QuantificationLevel
from apb_fasta.configuration import FastaAnnotationParameters
from cyclopts import App, Parameter
from loguru import logger
from pydantic import ValidationError

from apb_proteobench.api import (
    ScoredResult,
    VendorBenchmarkResult,
    benchmark_result,
    run_vendor_benchmark,
)
from apb_proteobench.integration import (
    ALL_ABUNDANCE_LAYERS,
    PRIMARY_LAYER,
    LayerSelection,
    NamedAbundanceLayer,
)
from apb_proteobench.io.result_performance import write_result_performance_bundle
from apb_proteobench.presentation import report_score

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
        Path | None,
        Parameter(help="ProteoBench module settings TOML"),
    ] = None
    output: Annotated[
        Path | None,
        Parameter(help="New scored APB2 .h5ad, .h5mu, .parquet, or .duckdb result"),
    ] = None
    software: Annotated[
        str | None,
        Parameter(help="Vendor software selecting APB2 parsing rules"),
    ] = None
    params_software: Annotated[
        str | None,
        Parameter(help="Software parser override for the parameter file"),
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
        str | None,
        Parameter(help="One named APB abundance layer to score"),
    ] = None
    x_only: Annotated[
        bool,
        Parameter(
            name="--x",
            negative=False,
            help="Score only the APB primary layer represented by X",
        ),
    ] = False
    result_performance: Annotated[
        Path | None,
        Parameter(help="Write result_performance.csv and sibling ProteoBot JSON"),
    ] = None
    strict: Annotated[
        bool,
        Parameter(negative=False, help="Promote APB2 layer-contract warnings to errors"),
    ] = False


DEFAULT_RUN_CLI_OPTIONS = RunCliOptions()


@app.command
def benchmark(
    source: Annotated[Path, Parameter(help="Existing APB2 result to annotate and score")],
    module: Annotated[Path, Parameter(help="ProteoBench module settings TOML")],
    target: Annotated[
        Path,
        Parameter(help="New scored APB2 result; format selected from its suffix"),
    ],
    /,
    *,
    layer: Annotated[
        str | None,
        Parameter(help="One named APB abundance layer to score"),
    ] = None,
    x_only: Annotated[
        bool,
        Parameter(
            name="--x",
            negative=False,
            help="Score only the APB primary layer represented by X",
        ),
    ] = False,
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
        selection = _layer_selection(
            layer,
            x_only,
            export_result_performance=result_performance is not None,
        )
    except _LayerSelectionUsageError as error:
        logger.error(str(error))
        return 2
    try:
        result = benchmark_result(source, module, target, selection=selection)
        _export_result_performance(result, result_performance)
    except (OSError, ValueError, ValidationError) as error:
        logger.error(str(error))
        return 1
    report_score(result, verbose=verbose)
    return 0


@app.command
def run(
    data: Annotated[Path, Parameter(help="Vendor result table or directory")],
    /,
    *fasta_paths: Annotated[
        Path,
        Parameter(help="One or more protein FASTA files"),
    ],
    options: Annotated[RunCliOptions, Parameter(name="*")] = DEFAULT_RUN_CLI_OPTIONS,
    verbose: Annotated[
        bool,
        Parameter(negative=False, help="Report detailed conversion and scoring diagnostics"),
    ] = False,
) -> int:
    """Write scored APB2 output and optional pMultiQC/ProteoBot export bundle.

    Quantitative aggregation is a separate step: run apb-aggregate between conversion
    and benchmarking when the scored level must be derived from a lower one. Use
    --level LEVEL to convert one quantification level. Scoring includes every declared
    abundance layer by default; use --x for only the APB primary/X layer.
    """
    if options.params is None:
        logger.error("pass --params PATH for the vendor search-parameter file")
        return 1
    if options.module is None:
        logger.error("pass --module PATH for the ProteoBench module settings")
        return 1
    if options.output is None:
        logger.error("pass --output PATH for the final APB2 result")
        return 1
    try:
        selection = _layer_selection(
            options.layer,
            options.x_only,
            export_result_performance=options.result_performance is not None,
        )
    except _LayerSelectionUsageError as error:
        logger.error(str(error))
        return 2
    try:
        result = run_vendor_benchmark(
            data,
            options.params,
            fasta_paths,
            options.module,
            options.output,
            level=options.level,
            software=options.software,
            parameters_software=options.params_software,
            checks="strict" if options.strict else "standard",
            selection=selection,
            fasta_parameters=FastaAnnotationParameters(
                protein_group_separator=options.protein_group_separator,
                matcher_backend=options.backend,
                il_equivalent=options.il_equivalent,
            ),
        )
        _export_result_performance(result.scored, options.result_performance)
    except (OSError, ValueError, ValidationError) as error:
        logger.error(str(error))
        return 1
    _report_vendor_benchmark(result, verbose=verbose)
    return 0


def _report_vendor_benchmark(result: VendorBenchmarkResult, /, *, verbose: bool) -> None:
    logger.info(
        "vendor={} software_version={}",
        result.software,
        result.software_version or "missing",
    )
    for level, coverage in result.fasta_reports.peptide_levels.items():
        logger.info(
            "level={} peptides_in_fasta={}/{} unmatched={}",
            level,
            coverage.matched_feature_count,
            coverage.feature_count,
            coverage.unmatched_feature_count,
        )
    report_score(result.scored, verbose=verbose)


class _LayerSelectionUsageError(ValueError):
    """The CLI received mutually exclusive layer-selection options."""


def _layer_selection(
    layer: str | None,
    x_only: bool,
    /,
    *,
    export_result_performance: bool = False,
) -> LayerSelection:
    if layer is not None and x_only:
        raise _LayerSelectionUsageError("--layer and --x are mutually exclusive")
    if export_result_performance and layer is None and not x_only:
        raise _LayerSelectionUsageError("--result-performance requires --x or --layer NAME")
    if layer is not None:
        return NamedAbundanceLayer(layer)
    if x_only:
        return PRIMARY_LAYER
    return ALL_ABUNDANCE_LAYERS


def _export_result_performance(result: ScoredResult, target: Path | None, /) -> None:
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
    datapoint = _proteobot_datapoint(result, selected.layer_name)
    written = write_result_performance_bundle(
        selected.analysis.diagnostics.legacy,
        datapoint,
        target,
    )
    logger.info("wrote pMultiQC input {}", written.csv)
    logger.info("wrote ProteoBot datapoint {}", written.proteobot_json)


def _proteobot_datapoint(result: ScoredResult, layer_name: str, /) -> dict[str, object]:
    selected = result.layers[layer_name]
    score_document: object = json.loads(selected.analysis.scores.model_dump_json())
    if not isinstance(score_document, dict):
        raise TypeError("ProteoBench scores did not serialize to a JSON object")
    parameters = result.search_parameters
    intermediate_hash = selected.analysis.scores.intermediate_hash
    software_name = _text_parameter(parameters, "software_name") or "unknown"
    datapoint: dict[str, object] = {
        "id": f"{software_name.replace(' ', '_')}_{intermediate_hash[:12]}",
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
    datapoint["submission_comments"] = (
        f"\n\nDataset URL: https://proteobench.cubimed.rub.de/datasets/{intermediate_hash}/"
    )
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


def main() -> int:
    """Console-script entry point."""
    result = app()
    return int(result) if result is not None else 0


if __name__ == "__main__":
    sys.exit(main())

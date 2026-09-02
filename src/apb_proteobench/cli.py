"""Cyclopts command line for vendor conversion, ProteoBench annotation, and scoring."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

from apb2.parserV2.vendor_parse_rules.schema.base import QuantificationLevel
from cyclopts import App, Parameter
from loguru import logger
from pydantic import ValidationError

from apb_proteobench.api import (
    ConvertedVendorResult,
    annotate_result,
    convert_vendor_result,
    score_result,
)
from apb_proteobench.presentation import report_score

app = App(
    name="apb-proteobench",
    help="Convert, annotate, and score ProteoBench experiments with APB2",
    help_on_error=True,
)

ANNDATA_SUFFIX = ".h5ad"
MUDATA_SUFFIX = ".h5mu"


@dataclass(frozen=True, slots=True)
class ConvertCliOptions:
    """Options for direct APB2 vendor conversion."""

    params: Path | None = None
    software: str | None = None
    params_software: str | None = None
    output: Path | None = None
    strict: bool = False


DEFAULT_CONVERT_CLI_OPTIONS = ConvertCliOptions()


@app.command
def convert(
    data: Path,
    level: QuantificationLevel | None = None,
    options: Annotated[ConvertCliOptions, Parameter(name="*")] = DEFAULT_CONVERT_CLI_OPTIONS,
) -> int:
    """Convert one vendor level to h5ad, or every compatible level to h5mu.

    --params supplies the required vendor search-parameter file. --software optionally selects
    and verifies the packaged vendor rules. --params-software independently selects the parameter
    parser. --output is a basename; this command appends .h5ad or .h5mu. --strict promotes layer
    contract warnings to errors.
    """
    if options.params is None:
        logger.error("pass --params PATH for the vendor search-parameter file")
        return 1
    output_suffix = MUDATA_SUFFIX if level is None else ANNDATA_SUFFIX
    if options.output is not None and options.output.suffix == output_suffix:
        logger.error(
            "--output must not already end in {}; apb-proteobench appends it, got {}",
            output_suffix,
            options.output,
        )
        return 2
    output = (
        data.with_suffix(output_suffix)
        if options.output is None
        else Path(f"{options.output}{output_suffix}")
    )
    try:
        result = convert_vendor_result(
            data,
            options.params,
            output,
            level=level,
            software=options.software,
            parameters_software=options.params_software,
            checks="strict" if options.strict else "standard",
        )
    except (OSError, ValueError, ValidationError) as error:
        logger.error(str(error))
        return 1
    _report_conversion(result)
    return 0


@app.command
def annotate(source: Path, module: Path, target: Path, /) -> int:
    """Annotate SOURCE using ProteoBench MODULE and write TARGET."""
    try:
        result = annotate_result(source, module, target)
    except (OSError, ValueError, ValidationError) as error:
        logger.error(str(error))
        return 1
    report = next(iter(result.reports.values()))
    logger.info(
        "annotated observations={} columns={} output={}",
        report.coverage.matched_observation_count,
        list(report.columns_added),
        target,
    )
    return 0


@app.command
def score(source: Path, target: Path, /, *, verbose: bool = False) -> int:
    """Score SOURCE from its embedded ProteoBench configuration and write TARGET."""
    try:
        result = score_result(source, target)
    except (OSError, ValueError, ValidationError) as error:
        logger.error(str(error))
        return 1
    report_score(result, verbose=verbose)
    return 0


def _report_conversion(result: ConvertedVendorResult) -> None:
    logger.info(
        "vendor={} software_version={}",
        result.software,
        result.software_version or "missing",
    )
    for level, parsed in result.parsed.levels.items():
        logger.info(
            "level={} shape=({}, {}) layers={}",
            level,
            parsed.obs.frame.height,
            parsed.var.frame.height,
            list(parsed.layers),
        )
    logger.info("wrote {}", result.output_path)


def main() -> int:
    """Console-script entry point."""
    result = app()
    return int(result) if result is not None else 0


if __name__ == "__main__":
    sys.exit(main())

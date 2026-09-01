"""Cyclopts command line for ProteoBench annotation and scoring."""

from __future__ import annotations

import sys
from pathlib import Path

from cyclopts import App
from loguru import logger
from pydantic import ValidationError

from apb_proteobench.api import annotate_result, score_result
from apb_proteobench.presentation import report_score

app = App(
    name="apb-proteobench",
    help="Annotate and score ProteoBench experiments stored as APB2 results",
    help_on_error=True,
)


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


def main() -> int:
    """Console-script entry point."""
    result = app()
    return int(result) if result is not None else 0


if __name__ == "__main__":
    sys.exit(main())

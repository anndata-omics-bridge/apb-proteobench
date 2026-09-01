"""Result-I/O workflows exposed to Python and the command line."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from apb2.annotation_extension import AnnotationResult
from apb2.result_facade import read_parsed_levels, write_parsed_levels

from apb_proteobench.annotation import ProteoBenchAnnotationParser
from apb_proteobench.configuration.schema import ModuleSettings
from apb_proteobench.integration import (
    ExtractedProteoBenchLevel,
    embedded_configuration,
    extract_level,
    persist_result,
)
from apb_proteobench.workflow import (
    DiagnosticMethod,
    MixedSpeciesDiagnostics,
    ProteoBenchCompatibleScoring,
    ProteoBenchResult,
    ScoringMethod,
    analyze_level,
)

_DEFAULT_DIAGNOSTICS = MixedSpeciesDiagnostics()
_DEFAULT_SCORING = ProteoBenchCompatibleScoring()


@dataclass(frozen=True, slots=True)
class ScoredResult:
    """Typed workflow evidence used by presenters and API clients."""

    input_path: Path
    output_path: Path
    configuration: ModuleSettings
    extracted: ExtractedProteoBenchLevel
    analysis: ProteoBenchResult


def annotate_result(source: Path, module: Path, target: Path, /) -> AnnotationResult:
    """Bind module samples/configuration to an APB result and persist a new result."""
    _require_new_target(source, target)
    parsed = read_parsed_levels(source)
    annotation = ProteoBenchAnnotationParser.from_path(module).parse(parsed)
    result = annotation.annotate()
    target.parent.mkdir(parents=True, exist_ok=True)
    write_parsed_levels(result.parsed, target)
    return result


def score_result(
    source: Path,
    target: Path,
    /,
    *,
    diagnostic_method: DiagnosticMethod = _DEFAULT_DIAGNOSTICS,
    scoring_method: ScoringMethod = _DEFAULT_SCORING,
) -> ScoredResult:
    """Score the configured APB level and persist a new result."""
    _require_new_target(source, target)
    parsed = read_parsed_levels(source)
    configuration = embedded_configuration(parsed)
    extracted = extract_level(parsed, configuration)
    analysis = analyze_level(
        extracted.calculation,
        configuration,
        diagnostic_method,
        scoring_method,
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    write_parsed_levels(persist_result(parsed, extracted, analysis), target)
    return ScoredResult(
        input_path=source,
        output_path=target,
        configuration=configuration,
        extracted=extracted,
        analysis=analysis,
    )


def _require_new_target(source: Path, target: Path) -> None:
    if source.resolve() == target.resolve():
        raise ValueError("output must differ from input")
    if target.exists():
        raise ValueError(f"output already exists: {target}")

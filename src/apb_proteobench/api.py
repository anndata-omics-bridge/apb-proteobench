"""Vendor conversion and result-I/O workflows exposed to Python and the CLI."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

from apb2.annotation_extension import AnnotationResult
from apb2.parserV2.compile import AnnDataOutput, ParseRuleCompiler, compile_mudata_parsers
from apb2.parserV2.detect_document import (
    DetectedRuleDocument,
    detect_rule_document,
    guess_software,
    search_parameter_evidence,
    software_slug,
)
from apb2.parserV2.parse_quant.data.parsed import (
    JsonValue,
    ParsedLevel,
    ParsedLevelName,
    ParsedLevels,
)
from apb2.parserV2.parse_quant.parameters.source import SingleFile
from apb2.parserV2.parse_rule_facade import PRODUCER, ParseRuleFacade
from apb2.parserV2.vendor_params.parsers.shared.model import Parameters
from apb2.parserV2.vendor_params.registry import parse_params
from apb2.parserV2.vendor_parse_rules.document import SearchParameterEvidence
from apb2.parserV2.vendor_parse_rules.schema.base import QuantificationLevel
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

type AnnDataChecks = Literal["standard", "strict"]


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
    extracted: ExtractedProteoBenchLevel
    analysis: ProteoBenchResult


def convert_vendor_result(
    data: Path,
    parameters_path: Path,
    target: Path,
    /,
    *,
    level: QuantificationLevel | None = None,
    software: str | None = None,
    parameters_software: str | None = None,
    checks: AnnDataChecks = "standard",
) -> ConvertedVendorResult:
    """Parse a vendor table with APB2's compiler and persist h5ad or h5mu.

    Args:
        data: Vendor result table.
        parameters_path: Vendor search-parameter file.
        target: Exact output path, ending in ``.h5ad`` for one level or ``.h5mu`` for all
            compatible levels.
        level: One quantification level. If omitted, parse every compatible level.
        software: Optional vendor slug used to select and verify the packaged rule document.
        parameters_software: Optional independent parameter-parser slug.
        checks: AnnData layer-contract validation level.

    Returns:
        The detected software metadata and parsed in-memory APB2 result.

    Raises:
        ValueError: The inputs cannot select or satisfy one packaged rule document, or the output
            suffix does not match the requested conversion.
    """
    _require_vendor_target(target, level)
    source, detected, parameters = _detect_vendor(
        data,
        parameters_path,
        software=software,
        parameters_software=parameters_software,
    )
    evidence = search_parameter_evidence(parameters)
    provenance = _vendor_provenance(detected, parameters, parameters_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    parsed = (
        _parse_all_levels(source, detected, evidence, provenance, target, checks=checks)
        if level is None
        else _parse_one_level(
            source,
            detected,
            evidence,
            provenance,
            level,
            target,
            checks=checks,
        )
    )
    return ConvertedVendorResult(
        input_path=data,
        parameters_path=parameters_path,
        output_path=target,
        software=detected.software,
        software_version=detected.version,
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


def _require_vendor_target(target: Path, level: QuantificationLevel | None) -> None:
    expected = ".h5mu" if level is None else ".h5ad"
    if target.suffix != expected:
        raise ValueError(f"output must end in {expected} for this conversion, got {target}")


def _detect_vendor(
    data: Path,
    parameters_path: Path,
    *,
    software: str | None,
    parameters_software: str | None,
) -> tuple[SingleFile, DetectedRuleDocument, Parameters]:
    source = SingleFile(path=data)
    requested_software = None if software is None else software_slug(software)
    parser_slug = parameters_software or requested_software or guess_software(source)
    if parser_slug is None:
        raise ValueError(f"could not auto-detect the vendor for {data}; pass --software SLUG")
    parameters = parse_params(parameters_path, software=software_slug(parser_slug))
    detected = detect_rule_document(parameters, source)
    if requested_software is not None and detected.software != requested_software:
        raise ValueError(
            f"software {requested_software!r} does not match the detected vendor "
            f"{detected.software!r}"
        )
    return source, detected, parameters


def _vendor_provenance(
    detected: DetectedRuleDocument,
    parameters: Parameters,
    parameters_path: Path,
) -> dict[str, JsonValue]:
    return {
        "rule_selection_method": (
            "software_version" if detected.version is not None else "columns"
        ),
        "search_parameters_version_status": (
            "missing" if parameters.software_version is None else "present"
        ),
        "search_parameters_path": str(parameters_path),
        "search_parameters": json.dumps(parameters.model_dump(mode="json")),
    }


def _parse_one_level(
    source: SingleFile,
    detected: DetectedRuleDocument,
    evidence: SearchParameterEvidence,
    provenance: dict[str, JsonValue],
    level: QuantificationLevel,
    target: Path,
    *,
    checks: AnnDataChecks,
) -> ParsedLevels:
    parser = ParseRuleCompiler(
        facade=ParseRuleFacade(detected.document, level, evidence),
        output=AnnDataOutput(checks=checks),
    ).compile(source)
    parsed_level = parser.parse()
    parsed_level.uns.update(provenance)
    parser.convert(parsed_level, target)
    return ParsedLevels(
        levels={level: parsed_level},
        uns={
            "produced_by": PRODUCER,
            **provenance,
            "quantification_levels": [level],
        },
    )


def _parse_all_levels(
    source: SingleFile,
    detected: DetectedRuleDocument,
    evidence: SearchParameterEvidence,
    provenance: dict[str, JsonValue],
    target: Path,
    *,
    checks: AnnDataChecks,
) -> ParsedLevels:
    parsers, writer = compile_mudata_parsers(
        document=detected.document,
        levels=detected.document.levels,
        parameter_evidence=evidence,
        source=source,
        checks=checks,
    )
    levels: dict[ParsedLevelName, ParsedLevel] = {}
    for parser in parsers:
        parsed_level = parser.parse()
        parsed_level.uns.update(provenance)
        levels[cast(ParsedLevelName, parser.level)] = parsed_level
    parsed = ParsedLevels(
        levels=levels,
        uns={
            "produced_by": PRODUCER,
            **provenance,
            "quantification_levels": list(levels),
        },
    )
    writer.write(parsed, target)
    return parsed

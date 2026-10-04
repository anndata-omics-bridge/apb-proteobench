"""Public in-memory ProteoBench analysis API."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from apb2.api import ParsedLevels

from apb_proteobench.annotation import ProteoBenchAnnotationParser
from apb_proteobench.configuration.entrapment import load_packaged_entrapment_module
from apb_proteobench.configuration.load import LoadedModule, load_module, load_packaged_module
from apb_proteobench.configuration.schema import ModuleSettings
from apb_proteobench.entrapment import EntrapmentAnalysisResult, EntrapmentAnalyzer
from apb_proteobench.integration import (
    ScoredLayerResult,
    SubmissionContent,
    diagnostics_slot,
    extract_layer,
    persist_results,
    select_layers,
    submission_content,
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


@dataclass(frozen=True, slots=True)
class ProteoBenchAnalysisResult:
    """A complete ProteoBench analysis and its scored APB2 result."""

    parsed: ParsedLevels
    configuration: ModuleSettings
    layers: dict[str, ScoredLayerResult]

    def submission(self, layer_name: str) -> SubmissionContent:
        """Return what one scored layer contributes to ProteoBench's submission files."""
        return submission_content(self.parsed, self.configuration, self.layers[layer_name])


class ProteoBenchAnalyzer:
    """Bind one validated ProteoBench module and its analysis collaborators."""

    __slots__ = (
        "_annotation_parser",
        "_configuration",
        "_diagnostic_method",
        "_layers",
        "_scoring_method",
    )

    def __init__(
        self,
        module: LoadedModule,
        layers: Sequence[str] | None = None,
        diagnostic_method: DiagnosticMethod = _DEFAULT_DIAGNOSTICS,
        scoring_method: ScoringMethod = _DEFAULT_SCORING,
    ) -> None:
        """Create a complete analyzer for one validated ProteoBench module.

        ``layers`` names the abundance layers to score; ``None`` scores every one.
        """
        self._annotation_parser = ProteoBenchAnnotationParser(module)
        self._configuration = module.settings
        self._layers = layers
        self._diagnostic_method = diagnostic_method
        self._scoring_method = scoring_method

    def analyze(self, parsed: ParsedLevels) -> ProteoBenchAnalysisResult:
        """Annotate and score canonical APB2 levels without physical I/O."""
        annotated = self._annotation_parser.parse(parsed).annotate()
        layers: dict[str, ScoredLayerResult] = {}
        for layer_name in select_layers(annotated, self._configuration, self._layers):
            selected = extract_layer(annotated, self._configuration, layer_name)
            analysis = analyze_level(
                selected.calculation,
                self._configuration,
                self._diagnostic_method,
                self._scoring_method,
            )
            layers[selected.layer_name] = ScoredLayerResult(
                level_name=selected.level_name,
                layer_name=selected.layer_name,
                diagnostics_slot=diagnostics_slot(selected.layer_name),
                roles=selected.roles,
                analysis=analysis,
            )
        return ProteoBenchAnalysisResult(
            parsed=persist_results(annotated, layers),
            configuration=self._configuration,
            layers=layers,
        )


__all__ = [
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

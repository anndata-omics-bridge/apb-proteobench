"""Public in-memory ProteoBench analysis API."""

from __future__ import annotations

from dataclasses import dataclass

from apb2.api import ParsedLevels

from apb_proteobench.annotation import ProteoBenchAnnotationParser
from apb_proteobench.configuration.load import LoadedModule
from apb_proteobench.configuration.schema import ModuleSettings
from apb_proteobench.integration import (
    ALL_ABUNDANCE_LAYERS,
    LayerSelection,
    ResolvedLayerSelection,
    ScoredLayerResult,
    diagnostics_slot,
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


@dataclass(frozen=True, slots=True)
class ProteoBenchAnalysisResult:
    """A complete ProteoBench analysis and its scored APB2 result."""

    parsed: ParsedLevels
    configuration: ModuleSettings
    selection: ResolvedLayerSelection
    layers: dict[str, ScoredLayerResult]


class ProteoBenchAnalyzer:
    """Bind one validated ProteoBench module and its analysis collaborators."""

    __slots__ = (
        "_annotation_parser",
        "_configuration",
        "_diagnostic_method",
        "_scoring_method",
        "_selection",
    )

    def __init__(
        self,
        module: LoadedModule,
        /,
        *,
        selection: LayerSelection = ALL_ABUNDANCE_LAYERS,
        diagnostic_method: DiagnosticMethod = _DEFAULT_DIAGNOSTICS,
        scoring_method: ScoringMethod = _DEFAULT_SCORING,
    ) -> None:
        """Create a complete analyzer for one validated ProteoBench module."""
        self._annotation_parser = ProteoBenchAnnotationParser(module)
        self._configuration = module.settings
        self._selection = selection
        self._diagnostic_method = diagnostic_method
        self._scoring_method = scoring_method

    def analyze(self, parsed: ParsedLevels, /) -> ProteoBenchAnalysisResult:
        """Annotate and score canonical APB2 levels without physical I/O."""
        annotated = self._annotation_parser.parse(parsed).annotate().parsed
        resolved = resolve_layer_selection(
            annotated,
            self._configuration,
            self._selection,
        )
        layers: dict[str, ScoredLayerResult] = {}
        for layer_name in resolved.layer_names:
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
            parsed=persist_results(annotated, resolved, layers),
            configuration=self._configuration,
            selection=resolved,
            layers=layers,
        )


__all__ = [
    "ProteoBenchAnalysisResult",
    "ProteoBenchAnalyzer",
]

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
    ProteinSource,
    ScoredLayerResult,
    SubmissionContent,
    diagnostics_slot,
    extract_layer,
    module_layer,
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
        "_proteins",
        "_scoring_method",
    )

    def __init__(
        self,
        module: LoadedModule,
        layers: Sequence[str] | None = None,
        proteins: ProteinSource = "matched",
        diagnostic_method: DiagnosticMethod = _DEFAULT_DIAGNOSTICS,
        scoring_method: ScoringMethod = _DEFAULT_SCORING,
    ) -> None:
        """Create a complete analyzer for one validated ProteoBench module.

        ``layers`` names the abundance layers to score; ``None`` scores every one, unless the
        module sets ``abundance_normalisation``: then ``None`` scores the layer apb-catalog's
        ``abundance`` set names for it, or the primary layer when the set names none.
        ``proteins`` picks which proteins decide species and contaminants: ``"matched"``, the
        FASTA proteins containing the peptide, or ``"vendor"``, the proteins the vendor reports,
        as ProteoBench decides. Both add the contaminants apb2 marked ``apb_Contaminant``.
        """
        self._annotation_parser = ProteoBenchAnnotationParser(module)
        self._configuration = module.settings
        self._layers = layers
        self._proteins: ProteinSource = proteins
        self._diagnostic_method = diagnostic_method
        self._scoring_method = scoring_method

    def analyze(self, parsed: ParsedLevels) -> ProteoBenchAnalysisResult:
        """Annotate and score canonical APB2 levels without physical I/O."""
        annotated = self._annotation_parser.parse(parsed).annotate()
        chosen = module_layer(annotated, self._configuration) if self._layers is None else None
        names = self._layers if chosen is None else (chosen.name,)
        layers: dict[str, ScoredLayerResult] = {}
        for layer_name in select_layers(annotated, self._configuration, names):
            selected = extract_layer(annotated, self._configuration, layer_name, self._proteins)
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
                proteins=self._proteins,
                analysis=analysis,
            )
        return ProteoBenchAnalysisResult(
            parsed=persist_results(annotated, layers, chosen),
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

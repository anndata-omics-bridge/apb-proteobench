"""Backend-independent ProteoBench calculation orchestration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from apb_proteobench.calculation.contracts import QuantitativeLevelInput
from apb_proteobench.calculation.intermediate import (
    IntermediateResult,
    align_runs,
    compute_intermediate,
)
from apb_proteobench.calculation.metrics import ProteoBenchScores, ScoreConfig, build_scores
from apb_proteobench.configuration.schema import ModuleSettings


class DiagnosticMethod(Protocol):
    """Calculate feature-aligned benchmark diagnostics."""

    def diagnose(
        self,
        level: QuantitativeLevelInput,
        configuration: ModuleSettings,
        /,
    ) -> IntermediateResult: ...

    def identity(self) -> dict[str, str]:
        """Return the persistable method identity."""
        ...


class ScoringMethod(Protocol):
    """Reduce one compatible diagnostic result to benchmark scores."""

    def score(
        self,
        diagnostics: IntermediateResult,
        configuration: ModuleSettings,
        /,
    ) -> ProteoBenchScores: ...

    def identity(self) -> dict[str, str]:
        """Return the persistable method identity."""
        ...


@dataclass(frozen=True, slots=True)
class MixedSpeciesDiagnostics:
    """HYE/HY diagnostics driven entirely by module species configuration."""

    def diagnose(
        self,
        level: QuantitativeLevelInput,
        configuration: ModuleSettings,
        /,
    ) -> IntermediateResult:
        design = align_runs(level.observations, configuration)
        return compute_intermediate(
            level.matrix,
            level.feature_ids,
            level.reported_proteins,
            configuration,
            design,
            level.level,
        )

    def identity(self) -> dict[str, str]:
        return {"name": "mixed-species-quantitative", "version": "1"}


@dataclass(frozen=True, slots=True)
class ProteoBenchCompatibleScoring:
    """ProteoBench 0.17-compatible aggregate scoring."""

    def score(
        self,
        diagnostics: IntermediateResult,
        configuration: ModuleSettings,
        /,
    ) -> ProteoBenchScores:
        general = configuration.general
        return build_scores(
            diagnostics.legacy,
            diagnostics.intermediate_hash,
            ScoreConfig(
                default_cutoff=general.default_cutoff_min_feature,
                max_nr_observed=general.max_nr_observed,
            ),
        )

    def identity(self) -> dict[str, str]:
        return {"name": "proteobench-compatible", "version": "1"}


@dataclass(frozen=True, slots=True)
class ProteoBenchResult:
    """Complete calculated result and the methods that produced it."""

    diagnostics: IntermediateResult
    scores: ProteoBenchScores
    diagnostic_method: dict[str, str]
    scoring_method: dict[str, str]


def analyze_level(
    level: QuantitativeLevelInput,
    configuration: ModuleSettings,
    diagnostic_method: DiagnosticMethod,
    scoring_method: ScoringMethod,
    /,
) -> ProteoBenchResult:
    """Run explicitly composed diagnostics and scoring for one level."""
    diagnostics = diagnostic_method.diagnose(level, configuration)
    return ProteoBenchResult(
        diagnostics=diagnostics,
        scores=scoring_method.score(diagnostics, configuration),
        diagnostic_method=diagnostic_method.identity(),
        scoring_method=scoring_method.identity(),
    )

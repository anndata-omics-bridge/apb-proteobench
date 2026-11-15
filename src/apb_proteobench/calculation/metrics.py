"""ProteoBench-compatible aggregate HYE/HY score metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict, Field, RootModel, model_validator
from scipy.stats import rankdata

PROTEOBENCH_COMPATIBILITY_VERSION = "0.18.0"
PROTEOBENCH_SOURCE_REVISION = "cd2a8f0034440bd337c9d34c41be0d36c08e19a7"

_HUMAN = "HUMAN"
_ABS_EPSILON = {"median": pl.col("abs_epsilon").median(), "mean": pl.col("abs_epsilon").mean()}

type ScoreMetric = int | float


class _ScoreModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, ser_json_inf_nan="null")


class ScoreConfig(_ScoreModel):
    """Module species and cutoff range used to build one ProteoBench score document."""

    species: tuple[str, ...] = Field(min_length=1)
    default_cutoff: int = Field(default=3, ge=1)
    max_nr_observed: int = Field(default=6, ge=1)

    @model_validator(mode="after")
    def _validate_cutoff(self) -> ScoreConfig:
        if self.default_cutoff > self.max_nr_observed:
            raise ValueError("default_cutoff must not exceed max_nr_observed")
        return self


class CutoffScores(RootModel[dict[str, ScoreMetric]]):
    """All compatible score metrics calculated at one observation cutoff."""

    model_config = ConfigDict(frozen=True, ser_json_inf_nan="null")


class ProteoBenchScores(_ScoreModel):
    """Typed storage contract for one ProteoBench score document."""

    results: dict[str, CutoffScores]
    proteobench_version: str
    median_abs_epsilon_global: float
    mean_abs_epsilon_global: float
    median_abs_epsilon_eq_species: float
    mean_abs_epsilon_eq_species: float
    median_abs_epsilon_precision_global: float
    mean_abs_epsilon_precision_global: float
    median_abs_epsilon_precision_eq_species: float
    mean_abs_epsilon_precision_eq_species: float
    nr_feature: int
    median_abs_epsilon_changed_global: float
    nr_quantified_changed: int
    dynamic_range_HUMAN: float
    median_abs_epsilon_HUMAN: float


def build_scores(
    intermediate: pd.DataFrame,
    config: ScoreConfig,
) -> ProteoBenchScores:
    """Compute the compatible score-only ProteoBench storage model."""
    inputs = _score_inputs(intermediate)
    unchanged = _unchanged_species(intermediate)
    cutoffs = {
        cutoff: inputs.filter(pl.col("nr_observed") >= cutoff)
        for cutoff in range(1, config.max_nr_observed + 1)
    }
    results = {
        str(cutoff): CutoffScores(
            _metrics_at_cutoff(intermediate, cutoff)
            | _changed_metrics(selected, config.species, unchanged)
            | _dynamic_ranges(selected, config.species)
        )
        for cutoff, selected in cutoffs.items()
    }
    selected = results[str(config.default_cutoff)]
    selected_metrics = selected.root
    return ProteoBenchScores(
        results=results,
        proteobench_version=PROTEOBENCH_COMPATIBILITY_VERSION,
        median_abs_epsilon_global=_required_float_metric(
            selected_metrics, "median_abs_epsilon_global"
        ),
        mean_abs_epsilon_global=_required_float_metric(selected_metrics, "mean_abs_epsilon_global"),
        median_abs_epsilon_eq_species=_required_float_metric(
            selected_metrics, "median_abs_epsilon_eq_species"
        ),
        mean_abs_epsilon_eq_species=_required_float_metric(
            selected_metrics, "mean_abs_epsilon_eq_species"
        ),
        median_abs_epsilon_precision_global=_required_float_metric(
            selected_metrics, "median_abs_epsilon_precision_global"
        ),
        mean_abs_epsilon_precision_global=_required_float_metric(
            selected_metrics, "mean_abs_epsilon_precision_global"
        ),
        median_abs_epsilon_precision_eq_species=_required_float_metric(
            selected_metrics, "median_abs_epsilon_precision_eq_species"
        ),
        mean_abs_epsilon_precision_eq_species=_required_float_metric(
            selected_metrics, "mean_abs_epsilon_precision_eq_species"
        ),
        nr_feature=_required_int_metric(selected_metrics, "nr_feature"),
        median_abs_epsilon_changed_global=_required_float_metric(
            selected_metrics, "median_abs_epsilon_changed_global"
        ),
        nr_quantified_changed=_required_int_metric(selected_metrics, "nr_quantified_changed"),
        dynamic_range_HUMAN=_species_metric(selected_metrics, f"dynamic_range_{_HUMAN}"),
        median_abs_epsilon_HUMAN=_species_metric(selected_metrics, f"median_abs_epsilon_{_HUMAN}"),
    )


def compute_roc_auc(frame: pd.DataFrame) -> float:
    """Compute binary ROC-AUC from absolute fold changes using average tie ranks."""
    required = {"species", "log2_A_vs_B", "log2_expectedRatio"}
    if frame.empty or not required <= set(frame.columns):
        return np.nan
    unchanged = _unchanged_species(frame)

    y_true = (frame["species"] != unchanged).to_numpy(dtype=np.int8)
    y_score = frame["log2_A_vs_B"].abs().to_numpy(dtype=np.float64)
    valid = ~np.isnan(y_score)
    y_true = y_true[valid]
    y_score = y_score[valid]
    if len(y_true) < 2 or len(np.unique(y_true)) < 2:
        return np.nan

    ranks = rankdata(y_score, method="average")
    positives = y_true == 1
    n_positive = int(np.count_nonzero(positives))
    n_negative = len(y_true) - n_positive
    rank_sum = float(ranks[positives].sum())
    return (rank_sum - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)


def _unchanged_species(frame: pd.DataFrame) -> object:
    """The species whose expected A:B ratio is closest to 1:1; ``None`` for an empty frame."""
    if frame.empty:
        return None
    species_ratios = frame[["species", "log2_expectedRatio"]].drop_duplicates()
    unchanged_index = species_ratios["log2_expectedRatio"].abs().idxmin()
    return species_ratios.loc[unchanged_index, "species"]


def _metrics_at_cutoff(frame: pd.DataFrame, cutoff: int) -> dict[str, ScoreMetric]:
    selected = frame[frame["nr_observed"] >= cutoff]
    return {
        **_species_metrics(selected),
        **_cv_metrics(selected),
        "variance_epsilon_global": float(selected["epsilon"].var()) if len(selected) else 0.0,
        "nr_feature": len(selected),
        "roc_auc": compute_roc_auc(selected),
    }


def _score_inputs(intermediate: pd.DataFrame) -> pl.DataFrame:
    """Select the changed-species and dynamic-range columns once; NaN becomes null."""
    columns = ["species", "nr_observed", "epsilon", "Intensity_mean_A", "Intensity_mean_B"]
    return pl.from_pandas(intermediate[columns], include_index=False).select(
        "species",
        "nr_observed",
        abs_epsilon=pl.col("epsilon").abs(),
        **{
            f"log10_{condition}": pl.col(f"Intensity_mean_{condition}")
            .clip(lower_bound=1e-10)
            .log10()
            for condition in "AB"
        },
    )


def _changed_metrics(
    selected: pl.DataFrame, species: tuple[str, ...], unchanged: object
) -> dict[str, ScoreMetric]:
    """Per-species counts and the error over the changed species, all but ``unchanged``.

    Ported from ProteoBench's plasma metrics, where the changed species are the spike-ins.
    As in ProteoBench, an empty row subset scores 0.0.
    """
    changed = selected if unchanged is None else selected.filter(pl.col("species") != unchanged)
    per_species = changed.group_by("species").agg(**_ABS_EPSILON)
    counts = dict(selected.group_by("species").len().iter_rows())
    metrics: dict[str, ScoreMetric] = {"nr_quantified_changed": changed.height}
    metrics |= {f"nr_quantified_{name}": counts.get(name, 0) for name in species}
    for stat, aggregate in _ABS_EPSILON.items():
        metrics[f"{stat}_abs_epsilon_changed_global"] = _subset_metric(changed, aggregate)
        metrics[f"{stat}_abs_epsilon_changed_eq_species"] = _subset_metric(
            per_species, pl.col(stat).mean()
        )
    return metrics


def _dynamic_ranges(selected: pl.DataFrame, species: tuple[str, ...]) -> dict[str, ScoreMetric]:
    """Each species' dynamic range per condition, A and B, and their mean.

    As in ProteoBench, an empty cutoff has zero dynamic range; a species without rows has NaN.
    """
    metrics: dict[str, ScoreMetric] = {}
    for name in species:
        rows = selected.filter(pl.col("species") == name)
        spread = [_dynamic_range(rows, condition) if selected.height else 0.0 for condition in "AB"]
        metrics[f"dynamic_range_{name}_A"], metrics[f"dynamic_range_{name}_B"] = spread
        metrics[f"dynamic_range_{name}"] = float(np.mean(spread))
    return metrics


def _subset_metric(frame: pl.DataFrame, aggregate: pl.Expr) -> float:
    return _as_float(frame.select(aggregate).item()) if frame.height else 0.0


def _dynamic_range(rows: pl.DataFrame, condition: str) -> float:
    """P90 minus P10 of log10 mean intensity, interpolated linearly as pandas does."""
    log10 = pl.col(f"log10_{condition}")
    return _as_float(
        rows.select(log10.quantile(0.9, "linear") - log10.quantile(0.1, "linear")).item()
    )


def _species_metrics(selected: pd.DataFrame) -> dict[str, float]:
    """Calculate global and equal-species metrics with one Polars grouping."""
    values = pl.from_pandas(
        selected[["species", "epsilon", "log2_A_vs_B"]], include_index=False
    ).with_columns(pl.col("epsilon", "log2_A_vs_B").fill_nan(None))
    values = values.with_columns(
        median_center=pl.col("log2_A_vs_B").median().over("species"),
        mean_center=pl.col("log2_A_vs_B").mean().over("species"),
    ).with_columns(
        median_precision=(pl.col("log2_A_vs_B") - pl.col("median_center")).abs(),
        mean_precision=(pl.col("log2_A_vs_B") - pl.col("mean_center")).abs(),
    )
    per_species = values.group_by("species").agg(
        pl.col("epsilon").abs().median().alias("median_abs_epsilon"),
        pl.col("epsilon").abs().mean().alias("mean_abs_epsilon"),
        pl.col("median_center").first(),
        pl.col("mean_center").first(),
        pl.col("median_precision").median(),
        pl.col("mean_precision").mean(),
    )
    global_metrics = values.select(
        pl.col("epsilon").abs().median().alias("median_abs_epsilon_global"),
        pl.col("epsilon").abs().mean().alias("mean_abs_epsilon_global"),
        pl.col("median_precision").median().alias("median_abs_epsilon_precision_global"),
        pl.col("mean_precision").mean().alias("mean_abs_epsilon_precision_global"),
    ).row(0, named=True)
    equal_species = per_species.select(
        pl.col("median_abs_epsilon").mean().alias("median_abs_epsilon_eq_species"),
        pl.col("mean_abs_epsilon").mean().alias("mean_abs_epsilon_eq_species"),
        pl.col("median_precision").mean().alias("median_abs_epsilon_precision_eq_species"),
        pl.col("mean_precision").mean().alias("mean_abs_epsilon_precision_eq_species"),
    ).row(0, named=True)
    metrics = {name: _as_float(value) for name, value in global_metrics.items()}
    metrics.update({name: _as_float(value) for name, value in equal_species.items()})
    for row in per_species.iter_rows(named=True):
        species = row["species"]
        if species is None:
            continue
        for name in ("median", "mean"):
            metrics[f"{name}_abs_epsilon_{species}"] = _as_float(row[f"{name}_abs_epsilon"])
            metrics[f"{name}_log2_empirical_{species}"] = _as_float(row[f"{name}_center"])
            metrics[f"{name}_abs_epsilon_precision_{species}"] = _as_float(row[f"{name}_precision"])
    return metrics


def _as_float(value: object) -> float:
    if value is None:
        return np.nan
    if not isinstance(value, (int, float)):
        raise TypeError(f"species metric is not numeric: {value!r}")
    return float(value)


def _cv_metrics(selected: pd.DataFrame) -> dict[str, float]:
    quantiles = selected[["CV_A", "CV_B"]].quantile([0.5, 0.75, 0.9, 0.95])
    averages = quantiles.mean(axis=1)
    return {
        "CV_median": float(averages.loc[0.50]),
        "CV_q75": float(averages.loc[0.75]),
        "CV_q90": float(averages.loc[0.90]),
        "CV_q95": float(averages.loc[0.95]),
    }


def _required_float_metric(metrics: dict[str, ScoreMetric], name: str) -> float:
    value = metrics[name]
    if isinstance(value, bool):
        raise TypeError(f"score metric {name!r} is not numeric")
    return float(value)


def _species_metric(metrics: dict[str, ScoreMetric], name: str) -> float:
    """A per-species metric; NaN when the species has no row at this cutoff."""
    return _required_float_metric(metrics, name) if name in metrics else np.nan


def _required_int_metric(metrics: dict[str, ScoreMetric], name: str) -> int:
    value = metrics[name]
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"score metric {name!r} is not an integer")
    return value

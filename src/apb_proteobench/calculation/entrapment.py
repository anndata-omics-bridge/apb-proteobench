"""ProteoBench-compatible entrapment labels and FDP metrics.

Ports ProteoBench's entrapment module (pair mapping in ``EntrapmentModule._apply_mapping``,
metrics in ``EntrapmentScores``) with one deliberate difference: precursors with equal
q-values share a rank. ProteoBench breaks such ties by vendor-file row order, which APB
results do not keep, so a tie never counts as an entrapment out-scoring its paired target.

Labels and pairs come from the entrapment FASTA instead of ProteoBench's pair file, which they
reproduce exactly: the FASTA lists each target directly before its entrapment, and the pair file
pairs each modified form with the partner form carrying the same modifications on the same
occurrence of each residue.
"""

from __future__ import annotations

import re
from typing import Literal

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict

type Category = Literal["valid", "invalid", "inconclusive"]

ENTRAPMENT_SOURCE_REVISION = "1147290715c32c454f31d2b20d7d13840a08c4e1"
"""The ProteoBench commit whose entrapment mapping and scores this module ports."""
TIE_RULE = "equal q-values share a rank"
PRECURSOR_COLUMNS = ("peptide", "sequence", "q_value")
_ALLOWED_MODIFICATIONS = frozenset({"carbamidomethyl", "oxidation"})
_MODIFICATION = re.compile(r"\[([^\[\]]*)\]")
_RESIDUE = re.compile(r"([A-Z])(?:\[([^\[\]]*)\])?")
_FIXED_THRESHOLDS = (0.001, 0.01, 0.05, 0.1, 1.0)
# 1 + 1 / entrapment fold, with one entrapment peptide per target (Wen et al. 2025, eq. 1).
_COMBINED_RATIO = 2.0
_EXAMPLES = 5


class EntrapmentError(ValueError):
    """The identified precursors cannot be labelled by the entrapment FASTA."""


class _ScoreModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, ser_json_inf_nan="null")


class ThresholdScores(_ScoreModel):
    """FDP estimates for the precursors at or below one q-value threshold."""

    lower_bound_FDP: float
    combined_FDP: float
    paired_FDP: float
    nr_id_features: int
    category_combined: Category
    category_paired: Category


class EntrapmentScores(_ScoreModel):
    """One ProteoBench entrapment score document."""

    nr_id_features: int
    reported_fdr_parsed_from_input: float
    combined_FDP: float
    lower_bound_FDP: float
    paired_FDP: float
    category_combined: Category
    category_paired: Category
    fdp_curve: dict[float, ThresholdScores]


def fasta_pairs(proteins: pl.DataFrame) -> pl.DataFrame:
    """Return the label and unmodified pair index of every entrapment-FASTA peptide.

    Args:
        proteins: The protein table of ProteoBench's entrapment FASTA, as
            ``FastaAnnotator.proteins`` holds it, in file order.

    Returns:
        One row per entry: ``peptide``, ``label`` (``target`` or ``entrapment``) and ``pair``.

    Raises:
        EntrapmentError: The frame lacks ``is_entrapment``, or two adjacent entries are not a
            target followed by its same-length entrapment.
    """
    if "is_entrapment" not in proteins.columns:
        raise EntrapmentError("the protein frame lacks protein_fasta's is_entrapment column")
    pairs = proteins.select(
        peptide="sequence",
        label=pl.when("is_entrapment").then(pl.lit("entrapment")).otherwise(pl.lit("target")),
        pair=pl.int_range(pl.len(), dtype=pl.Int64) // 2,
    )
    broken = (
        pairs.group_by("pair")
        .agg(
            labels=pl.col("label").str.join(","),
            lengths=pl.col("peptide").str.len_chars().n_unique(),
        )
        .filter((pl.col("labels") != "target,entrapment") | (pl.col("lengths") != 1))
    )
    if broken.height:
        raise EntrapmentError(
            f"{broken.height} adjacent FASTA entry pairs are not a target followed by its "
            f"same-length entrapment; is this ProteoBench's entrapment FASTA? First pair: "
            f"{broken.get_column('pair').min()}"
        )
    return pairs


def label_precursors(
    precursors: pl.DataFrame,
    pairs: pl.DataFrame,
    *,
    max_missing_fraction: float = 0.01,
) -> pl.DataFrame:
    """Drop precursors of unknown peptides and label the rest as target or entrapment.

    Args:
        precursors: One row per precursor: ``peptide`` (stripped), ``sequence`` (modified,
            with ProteoBench modification names such as ``M[Oxidation]``) and ``q_value``.
        pairs: The pair table from :func:`fasta_pairs`.
        max_missing_fraction: Largest tolerated fraction of peptides absent from ``pairs``.

    Returns:
        The kept precursors with ``label`` and ``pair``: the unmodified pair index followed by
        each modification and the occurrence of the residue it sits on.

    Raises:
        EntrapmentError: Too many peptides are unknown, or a sequence carries a modification
            other than oxidation or carbamidomethylation.
    """
    peptides = precursors.get_column("peptide").unique()
    missing = peptides.filter(~peptides.is_in(pairs.get_column("peptide").implode()))
    if missing.len() > max_missing_fraction * peptides.len():
        raise EntrapmentError(
            f"{missing.len()} of {peptides.len()} identified peptides are absent from the "
            f"entrapment FASTA (limit {max_missing_fraction:.0%}). Search the pre-digested "
            f"entrapment FASTA without enzymatic cleavage. First: {_first(missing)}"
        )
    kept = precursors.filter(~pl.col("peptide").is_in(missing.implode()))
    unsupported = [
        sequence
        for sequence in kept.get_column("sequence").unique().sort()
        if any(
            mod.strip().lower() not in _ALLOWED_MODIFICATIONS
            for mod in _MODIFICATION.findall(sequence)
        )
    ]
    if unsupported:
        raise EntrapmentError(
            f"{len(unsupported)} sequence(s) carry a modification other than Oxidation or "
            f"Carbamidomethyl, which ProteoBench's pairs do not cover. First: "
            f"{', '.join(unsupported[:_EXAMPLES])}"
        )
    sites = {sequence: _sites(sequence) for sequence in kept.get_column("sequence").unique()}
    return kept.join(pairs, on="peptide", how="left", maintain_order="left").with_columns(
        pl.format("{}{}", "pair", pl.col("sequence").replace_strict(sites)).alias("pair")
    )


def _sites(sequence: str) -> str:
    """Name each modification by residue and occurrence, e.g. ``:M2[oxidation]``."""
    seen: dict[str, int] = {}
    sites: list[str] = []
    for residue, modification in _RESIDUE.findall(sequence):
        seen[residue] = seen.get(residue, 0) + 1
        if modification:
            sites.append(f":{residue}{seen[residue]}[{modification.strip().lower()}]")
    return "".join(sorted(sites))


def score_entrapment(labelled: pl.DataFrame, *, intervals: int = 10) -> EntrapmentScores:
    """Compute ProteoBench's entrapment metrics for labelled precursors.

    Args:
        labelled: The output of :func:`label_precursors`, one row per precursor.
        intervals: Number of evenly spaced q-value thresholds on the FDP curve.

    Returns:
        The lower-bound, combined and paired FDP, their categories against the largest
        reported q-value, and the same estimates at each curve threshold.

    Raises:
        EntrapmentError: No precursor is labelled, or a q-value is missing.
    """
    if labelled.height == 0:
        raise EntrapmentError("no precursor carries an entrapment label")
    if labelled.get_column("q_value").null_count():
        raise EntrapmentError("every precursor needs a q-value")
    ranked = labelled.with_columns(pl.col("q_value").rank(method="min").alias("score"))
    reported_fdr = float(ranked.select(pl.col("q_value").max()).item())
    whole = _threshold_scores(ranked, reported_fdr)
    return EntrapmentScores(
        nr_id_features=whole.nr_id_features,
        reported_fdr_parsed_from_input=reported_fdr,
        combined_FDP=whole.combined_FDP,
        lower_bound_FDP=whole.lower_bound_FDP,
        paired_FDP=whole.paired_FDP,
        category_combined=whole.category_combined,
        category_paired=whole.category_paired,
        fdp_curve=_fdp_curve(ranked, reported_fdr, intervals),
    )


def _fdp_curve(ranked: pl.DataFrame, max_q: float, intervals: int) -> dict[float, ThresholdScores]:
    fixed = [threshold for threshold in _FIXED_THRESHOLDS if threshold <= max_q]
    spaced = np.linspace(max_q / intervals, max_q, intervals).tolist()
    curve: dict[float, ThresholdScores] = {}
    for threshold in sorted(set(fixed + spaced)):
        subset = ranked.filter(pl.col("q_value") <= threshold)
        if subset.filter(pl.col("label") == "target").height:
            curve[round(float(threshold), 8)] = _threshold_scores(subset, threshold)
    return curve


def _threshold_scores(ranked: pl.DataFrame, fdr: float) -> ThresholdScores:
    entrapments = ranked.filter(pl.col("label") == "entrapment")
    targets = ranked.filter(pl.col("label") == "target")
    identified = entrapments.height + targets.height
    lower = entrapments.height / identified
    combined = entrapments.height * _COMBINED_RATIO / identified
    paired = _paired_fdp(entrapments, targets) / identified
    return ThresholdScores(
        lower_bound_FDP=lower,
        combined_FDP=combined,
        paired_FDP=paired,
        nr_id_features=ranked.height,
        category_combined=_category(lower, combined, fdr),
        category_paired=_category(lower, paired, fdr),
    )


def _paired_fdp(entrapments: pl.DataFrame, targets: pl.DataFrame) -> int:
    """Return the paired-FDP numerator of Wen et al. 2025, eq. 2.

    Every identified entrapment counts once, again when its paired target is not
    identified, and twice more when it ranks strictly better than that target.
    """
    best_entrapment = entrapments.group_by("pair").agg(pl.col("score").max().alias("entrapment"))
    best_target = targets.group_by("pair").agg(pl.col("score").max().alias("target"))
    pairs = best_entrapment.join(best_target, on="pair", how="left")
    without_target = pairs.get_column("target").null_count()
    better = pairs.filter(pl.col("entrapment") < pl.col("target")).height
    return entrapments.height + without_target + 2 * better


def _category(lower: float, upper: float, fdr: float) -> Category:
    if upper <= fdr:
        return "valid"
    if lower > fdr:
        return "invalid"
    return "inconclusive"


def _first(values: pl.Series) -> str:
    return ", ".join(values.sort().head(_EXAMPLES).to_list())

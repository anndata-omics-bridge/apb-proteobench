"""Entrapment labels and FDP metrics against hand-computed values."""

from __future__ import annotations

import polars as pl
import pytest

from apb_proteobench.calculation.entrapment import (
    EntrapmentError,
    fasta_pairs,
    label_precursors,
    score_entrapment,
)

# The entrapment FASTA in file order: each target directly before its entrapment.
PROTEINS = pl.DataFrame(
    {
        "sequence": ["AAAK", "CCCK", "DDDK", "EEEK", "FFMK", "GGMK"],
        "is_entrapment": [False, True, False, True, False, True],
    }
)


@pytest.fixture
def pairs() -> pl.DataFrame:
    return fasta_pairs(PROTEINS)


def _precursors(*rows: tuple[str, str, float]) -> pl.DataFrame:
    return pl.DataFrame(rows, schema=["peptide", "sequence", "q_value"], orient="row")


# Pair 0: entrapment worse than its target; pair 1: tied; pair 2: entrapment better.
SIX = _precursors(
    ("AAAK", "AAAK", 0.001),
    ("CCCK", "CCCK", 0.002),
    ("DDDK", "DDDK", 0.003),
    ("EEEK", "EEEK", 0.003),
    ("FFMK", "FFM[Oxidation]K", 0.004),
    ("GGMK", "GGM[Oxidation]K", 0.0005),
)


def test_fasta_neighbours_are_pairs_and_entrapments_are_named(pairs: pl.DataFrame) -> None:
    assert pairs.rows()[:2] == [("AAAK", "target", 0), ("CCCK", "entrapment", 0)]
    assert pairs.get_column("pair").to_list() == [0, 0, 1, 1, 2, 2]


def test_a_fasta_that_does_not_alternate_target_and_entrapment_is_refused() -> None:
    swapped = PROTEINS.with_columns(pl.col("is_entrapment").shift(1, fill_value=True))
    with pytest.raises(EntrapmentError, match="3 adjacent FASTA entry pairs"):
        fasta_pairs(swapped)
    with pytest.raises(EntrapmentError, match="lacks protein_fasta's is_entrapment"):
        fasta_pairs(PROTEINS.drop("is_entrapment"))


def test_modified_forms_pair_by_modification_and_residue_occurrence(pairs: pl.DataFrame) -> None:
    labelled = label_precursors(
        _precursors(
            ("FFMK", "FFMK", 0.001),
            ("FFMK", "FFM[Oxidation]K", 0.001),
            ("GGMK", "GGM[Oxidation]K", 0.001),
        ),
        pairs,
    )

    assert labelled.get_column("pair").to_list() == ["2", "2:M1[oxidation]", "2:M1[oxidation]"]


def test_scores_match_hand_computed_values(pairs: pl.DataFrame) -> None:
    scores = score_entrapment(label_precursors(SIX, pairs))

    assert scores.nr_id_features == 6
    assert scores.reported_fdr_parsed_from_input == 0.004
    assert scores.lower_bound_FDP == 3 / 6
    assert scores.combined_FDP == 6 / 6
    # Three entrapments, none without its target, one ranked better: (3 + 0 + 2) / 6.
    assert scores.paired_FDP == 5 / 6
    assert (scores.category_combined, scores.category_paired) == ("invalid", "invalid")
    assert list(scores.fdp_curve) == [
        0.001,
        0.0012,
        0.0016,
        0.002,
        0.0024,
        0.0028,
        0.0032,
        0.0036,
        0.004,
    ], "thresholds without an identified target are left out"
    first = scores.fdp_curve[0.001]
    # AAAK and GGM[Oxidation]K, whose target FFM[Oxidation]K is not identified yet:
    # (1 + 1 + 0) / 2.
    assert (first.nr_id_features, first.lower_bound_FDP, first.paired_FDP) == (2, 0.5, 1.0)
    assert scores.fdp_curve[0.004].paired_FDP == scores.paired_FDP


def test_a_tie_is_not_an_entrapment_out_scoring_its_target(pairs: pl.DataFrame) -> None:
    tied = _precursors(("EEEK", "EEEK", 0.01), ("DDDK", "DDDK", 0.01))

    scores = score_entrapment(label_precursors(tied, pairs))

    assert scores.paired_FDP == (1 + 0 + 0) / 2


def test_too_many_unknown_peptides_are_refused(pairs: pl.DataFrame) -> None:
    with_unknown = pl.concat([SIX, _precursors(("ZZZK", "ZZZK", 0.001))])

    with pytest.raises(EntrapmentError, match="1 of 7 identified peptides are absent"):
        label_precursors(with_unknown, pairs)
    tolerated = label_precursors(with_unknown, pairs, max_missing_fraction=0.5)
    assert tolerated.get_column("peptide").to_list() == SIX.get_column("peptide").to_list()


def test_other_modifications_are_refused(pairs: pl.DataFrame) -> None:
    with pytest.raises(EntrapmentError, match="other than Oxidation"):
        label_precursors(_precursors(("AAAK", "AAAK[Phospho]", 0.001)), pairs)


def test_scoring_needs_labelled_precursors_with_q_values(pairs: pl.DataFrame) -> None:
    with pytest.raises(EntrapmentError, match="no precursor"):
        score_entrapment(label_precursors(SIX.clear(), pairs))
    missing_q = SIX.with_columns(pl.lit(None, dtype=pl.Float64).alias("q_value"))
    with pytest.raises(EntrapmentError, match="needs a q-value"):
        score_entrapment(label_precursors(missing_q, pairs))

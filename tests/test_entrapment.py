"""Entrapment labels and FDP metrics against hand-computed values."""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from apb_proteobench.calculation.entrapment import (
    EntrapmentError,
    label_precursors,
    read_pairs,
    score_entrapment,
)

PAIRS = """sequence\tdecoy\tproteins\tpeptide_type\tpeptide_pair_index
AAAK\tNo\tsp|P1|ONE_HUMAN\ttarget\t0
CCCK\tNo\tsp|P1_p_target|ONE_HUMAN_p_target\tp_target\t0
DDDK\tNo\tsp|P2|TWO_HUMAN\ttarget\t1
EEEK\tNo\tsp|P2_p_target|TWO_HUMAN_p_target\tp_target\t1
FFMK\tNo\tsp|P3|THREE_HUMAN\ttarget\t2
FFM[Oxidation]K\tNo\tsp|P3|THREE_HUMAN\ttarget\t2
GGMK\tNo\tsp|P3_p_target|THREE_HUMAN_p_target\tp_target\t2
AAAK\tNo\tsp|P9|NINE_HUMAN\ttarget\t9
"""


@pytest.fixture
def pairs(tmp_path: Path) -> pl.DataFrame:
    path = tmp_path / "pairs.txt.gz"
    path.write_bytes(gzip.compress(PAIRS.encode()))
    return read_pairs(path)


def _precursors(*rows: tuple[str, str, float]) -> pl.DataFrame:
    return pl.DataFrame(rows, schema=["peptide", "sequence", "q_value"], orient="row")


# Pair 0: entrapment worse than its target; pair 1: tied; pair 2: entrapment better.
SIX = _precursors(
    ("AAAK", "AAAK", 0.001),
    ("CCCK", "CCCK", 0.002),
    ("DDDK", "DDDK", 0.003),
    ("EEEK", "EEEK", 0.003),
    ("FFMK", "FFM[Oxidation]K", 0.004),
    ("GGMK", "GGMK", 0.0005),
)


def test_pairs_keep_the_first_row_per_sequence_and_name_entrapments(pairs: pl.DataFrame) -> None:
    assert pairs.row(0) == ("AAAK", "target", 0)
    assert pairs.filter(pl.col("label") == "entrapment").get_column("sequence").to_list() == [
        "CCCK",
        "EEEK",
        "GGMK",
    ]
    assert pairs.height == 7


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
    # AAAK and GGMK: GGMK's target FFMK is not identified yet: (1 + 1 + 0) / 2.
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


def test_other_modifications_and_unknown_modified_forms_are_refused(pairs: pl.DataFrame) -> None:
    with pytest.raises(EntrapmentError, match="other than Oxidation"):
        label_precursors(_precursors(("AAAK", "AAAK[Phospho]", 0.001)), pairs)
    with pytest.raises(EntrapmentError, match="absent from the pair file although"):
        label_precursors(_precursors(("GGMK", "GGM[Oxidation]K", 0.001)), pairs)


def test_scoring_needs_labelled_precursors_with_q_values(pairs: pl.DataFrame) -> None:
    with pytest.raises(EntrapmentError, match="no precursor"):
        score_entrapment(label_precursors(SIX.clear(), pairs))
    missing_q = SIX.with_columns(pl.lit(None, dtype=pl.Float64).alias("q_value"))
    with pytest.raises(EntrapmentError, match="needs a q-value"):
        score_entrapment(label_precursors(missing_q, pairs))

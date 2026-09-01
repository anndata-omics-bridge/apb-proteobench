"""Storage-independent ProteoBench calculation compatibility."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from apb_proteobench.calculation.contracts import QuantitativeLevelInput
from apb_proteobench.calculation.intermediate import align_runs, compute_intermediate
from apb_proteobench.calculation.mapping import render_proteobench_features
from apb_proteobench.calculation.metrics import ScoreConfig, build_scores, compute_roc_auc
from apb_proteobench.configuration.load import load_module
from apb_proteobench.workflow import (
    MixedSpeciesDiagnostics,
    ProteoBenchCompatibleScoring,
    analyze_level,
)
from conftest import matrix_values, module_settings, quantitative_input

GOLDEN = Path(__file__).parent / "data" / "small_legacy_intermediate.txt"
GOLDEN_HASH = "9077847f733c12b1297a4928a0b4c509e50e4ed9"


def test_hye_intermediate_matches_legacy_golden() -> None:
    inputs = quantitative_input()
    result = MixedSpeciesDiagnostics().diagnose(inputs, module_settings())

    pd.testing.assert_frame_equal(
        result.legacy,
        pd.read_csv(GOLDEN, index_col=0),
        check_dtype=False,
    )
    assert result.intermediate_hash == GOLDEN_HASH
    assert result.varm["included"].tolist() == [True, True, True, False, False, False]
    assert result.protein_mapping.accession_mapper.entries == 38_233
    assert result.protein_mapping.accession_mapper.sha256 == (
        "032034e2f9bea3fc41290c7461417280b1d37cec41ee8ef9a44c250781a4b997"
    )


def test_hy_uses_the_same_configuration_driven_calculation() -> None:
    result = analyze_level(
        quantitative_input(hy=True),
        module_settings(hy=True),
        MixedSpeciesDiagnostics(),
        ProteoBenchCompatibleScoring(),
    )

    assert set(result.diagnostics.legacy["species"]) == {"HUMAN", "YEAST"}
    assert result.scores.nr_feature == 3


def test_single_cell_hy_module_matches_hand_computed_ratios() -> None:
    configuration = load_module(Path(__file__).parent / "data" / "dia_singlecell.toml").settings
    samples = configuration.samples
    inputs = QuantitativeLevelInput(
        observations=pd.DataFrame(
            {
                "sample_name": [sample.sample_name for sample in samples],
                "condition": [sample.condition for sample in samples],
            }
        ),
        matrix=np.asarray([[12.0, 2.0]] * 3 + [[10.0, 10.0]] * 3),
        feature_ids=pd.Index(["H/2", "Y/2"]),
        reported_proteins=pd.Series(["P1_HUMAN", "P2_YEAST"]),
        level="ion",
    )

    result = analyze_level(
        inputs,
        configuration,
        MixedSpeciesDiagnostics(),
        ProteoBenchCompatibleScoring(),
    )

    assert set(result.diagnostics.legacy["species"]) == {"HUMAN", "YEAST"}
    assert result.diagnostics.legacy["epsilon"].abs().max() == pytest.approx(0.0, abs=2e-7)
    assert result.scores.nr_feature == 2
    assert result.scores.results["1"].root["roc_auc"] == 1.0


def test_dense_and_sparse_diagnostics_are_equal() -> None:
    inputs = quantitative_input()
    configuration = module_settings()
    design = align_runs(inputs.observations, configuration)
    dense = compute_intermediate(
        inputs.matrix,
        inputs.feature_ids,
        inputs.reported_proteins,
        configuration,
        design,
        "ion",
    )
    sparse_result = compute_intermediate(
        sparse.csr_matrix(np.nan_to_num(matrix_values(), nan=0.0)),
        inputs.feature_ids,
        inputs.reported_proteins,
        configuration,
        design,
        "ion",
    )
    pd.testing.assert_frame_equal(dense.varm, sparse_result.varm)
    pd.testing.assert_frame_equal(dense.legacy, sparse_result.legacy)


def test_score_names_cutoffs_and_roc_edge_cases() -> None:
    result = MixedSpeciesDiagnostics().diagnose(quantitative_input(), module_settings())
    scores = build_scores(result.legacy, result.intermediate_hash, ScoreConfig())
    tied = pd.DataFrame(
        {
            "species": ["HUMAN", "HUMAN", "YEAST", "YEAST"],
            "log2_A_vs_B": [0.0, 1.0, 0.0, 1.0],
            "log2_expectedRatio": [0.0, 0.0, 1.0, 1.0],
        }
    )

    assert list(scores.results) == ["1", "2", "3", "4", "5", "6"]
    assert scores.results["1"].root["nr_feature"] == 3
    assert scores.results["1"].root["roc_auc"] == 1.0
    assert np.isnan(scores.results["5"].root["CV_median"])
    assert compute_roc_auc(tied) == 0.5
    assert np.isnan(compute_roc_auc(tied[tied["species"] == "HUMAN"]))


def test_feature_rendering_uses_apb_canonical_unimod_vocabulary() -> None:
    features = pd.Series(["AC[UNIMOD:4]DM[UNIMOD:35]/2", "PEPTIDEC-[UNIMOD:4]/2"])

    assert render_proteobench_features(features).tolist() == [
        "AC[Carbamidomethyl]DM[Oxidation]/2",
        "PEPTIDEC-[Carbamidomethyl]/2",
    ]
    assert render_proteobench_features(
        features, drop_final_residue_modifications=True
    ).tolist() == ["AC[Carbamidomethyl]DM/2", "PEPTIDEC/2"]


def test_alignment_rejects_conflicting_design() -> None:
    inputs = quantitative_input()
    inputs.observations.loc[0, "condition"] = "B"

    with pytest.raises(ValueError, match="does not match module condition"):
        MixedSpeciesDiagnostics().diagnose(inputs, module_settings())

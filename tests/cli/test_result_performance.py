"""ProteoBench result-performance CSV output."""

from __future__ import annotations

import errno
import json
import os
from dataclasses import replace
from pathlib import Path
from typing import TextIO

import numpy as np
import pandas as pd
import pytest

from apb_proteobench.cli.result_performance import (
    content_hash,
    write_result_performance,
    write_result_performance_bundle,
)
from apb_proteobench.integration import SubmissionContent

GOLDEN = Path(__file__).parents[1] / "data" / "small_legacy_intermediate.txt"


def _ion_intermediate() -> pd.DataFrame:
    return pd.read_csv(GOLDEN, index_col=0)


def _submission_content() -> SubmissionContent:
    return SubmissionContent(
        matrix=np.array([[1.0, np.nan], [2.0, -0.0]]),
        feature_ids=pd.Index(["B/2", "A/2"]),
        reported_proteins=pd.Series(["P2_YEAST", "P1_HUMAN"]),
        raw_files=("run_b", "run_a"),
        conditions=("B", "A"),
        settings={"species_mapper": {"_HUMAN": "HUMAN", "_YEAST": "YEAST"}},
    )


def test_public_content_hash_is_stable() -> None:
    assert content_hash(_submission_content()) == content_hash(_submission_content())
    assert len(content_hash(_submission_content())) == 64


def _proteobot_datapoint() -> dict[str, object]:
    return {
        "software_name": "Synthetic",
        "software_version": "1.0",
        "search_engine": None,
        "search_engine_version": None,
        "ident_fdr_psm": 0.01,
        "ident_fdr_peptide": None,
        "ident_fdr_protein": 0.01,
        "enable_match_between_runs": False,
        "precursor_mass_tolerance": "[-10 ppm, 10 ppm]",
        "fragment_mass_tolerance": None,
        "enzyme": "Trypsin",
        "allowed_miscleavages": 2,
        "min_peptide_length": 7,
        "max_peptide_length": None,
        "is_temporary": False,
        "results": {"1": {"nr_feature": 3}},
        "median_abs_epsilon_global": 0.1,
        "mean_abs_epsilon_global": 0.2,
        "median_abs_epsilon_eq_species": 0.3,
        "mean_abs_epsilon_eq_species": 0.4,
        "median_abs_epsilon_precision_global": 0.5,
        "mean_abs_epsilon_precision_global": 0.6,
        "median_abs_epsilon_precision_eq_species": 0.7,
        "mean_abs_epsilon_precision_eq_species": 0.8,
        "nr_feature": 3,
        "comments": "",
        "proteobench_version": "0.17.0",
        "old_new": "new",
        "semi_enzymatic": False,
        "fixed_mods": "C[Carbamidomethyl]",
        "variable_mods": None,
        "max_mods": 3,
        "min_precursor_charge": 2,
        "max_precursor_charge": 6,
        "quantification_method": None,
        "protein_inference": None,
        "abundance_normalization_ions": None,
        "postprocessing_performed": False,
        "postprocessing_description": None,
        "submission_comments": "",
    }


def test_write_result_performance_round_trips_legacy_golden(tmp_path: Path) -> None:
    expected = _ion_intermediate()
    target = tmp_path / "pmultiqc" / "result_performance.csv"

    written = write_result_performance(expected, target)

    actual = pd.read_csv(written)
    pd.testing.assert_frame_equal(actual, expected.reset_index(drop=True), check_dtype=False)
    assert list(actual.columns) == list(expected.columns)
    assert not target.read_text(encoding="utf-8").splitlines()[0].startswith(",")
    assert list(target.parent.glob(".*.tmp")) == []


def test_write_result_performance_refuses_wrong_name_and_non_ion_table(
    tmp_path: Path,
) -> None:
    frame = _ion_intermediate()

    with pytest.raises(ValueError, match="must be named"):
        write_result_performance(frame, tmp_path / "other.csv")
    with pytest.raises(ValueError, match="ion-level"):
        write_result_performance(
            frame.rename(columns={"precursor ion": "peptidoform"}),
            tmp_path / "result_performance.csv",
        )


def test_write_result_performance_refuses_overwrite(tmp_path: Path) -> None:
    target = tmp_path / "result_performance.csv"
    target.write_text("owned by user\n", encoding="utf-8")

    with pytest.raises(ValueError, match="already exists"):
        write_result_performance(_ion_intermediate(), target)

    assert target.read_text(encoding="utf-8") == "owned by user\n"


def test_write_result_performance_removes_failed_temporary_file(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    target = tmp_path / "result_performance.csv"

    def fail_after_partial_write(
        _frame: pd.DataFrame,
        handle: TextIO,
        **_kwargs: object,
    ) -> None:
        handle.write("partial")
        raise OSError("write failed")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_after_partial_write)

    with pytest.raises(OSError, match="write failed"):
        write_result_performance(_ion_intermediate(), target)

    assert not target.exists()
    assert list(tmp_path.glob(".*.tmp")) == []


def test_write_result_performance_bundle_writes_proteobot_hash_json(tmp_path: Path) -> None:
    target = tmp_path / "bundle" / "result_performance.csv"
    datapoint = _proteobot_datapoint()
    content = _submission_content()

    written = write_result_performance_bundle(
        _ion_intermediate(), datapoint, target, content=content
    )

    assert written.csv == target
    digest = written.proteobot_json.stem
    assert written.proteobot_json == target.with_name(f"{digest}.json")
    assert digest == "3c52330156ef084849e0a7993cb986a76fc4ac6b2c9631f88ecdab93e75106ad"
    assert json.loads(written.proteobot_json.read_text(encoding="utf-8")) == {
        **datapoint,
        "id": f"Synthetic_{digest[:12]}",
        "intermediate_hash": digest,
        "submission_comments": (
            f"\n\nDataset URL: https://proteobench.cubimed.rub.de/datasets/{digest}/"
        ),
    }
    assert list(target.parent.glob(".*.tmp")) == []


def test_content_hash_follows_values_and_identities_not_axis_order(tmp_path: Path) -> None:
    content = _submission_content()
    reordered = SubmissionContent(
        matrix=content.matrix[[1, 0]][:, [1, 0]],
        feature_ids=pd.Index(list(reversed(content.feature_ids.tolist()))),
        reported_proteins=content.reported_proteins.iloc[np.array([1, 0])].reset_index(drop=True),
        raw_files=("run_a", "run_b"),
        conditions=("A", "B"),
        settings=content.settings,
    )
    changed_values = SubmissionContent(
        matrix=np.array([[1.0, np.nan], [2.0, 0.5]]),
        feature_ids=content.feature_ids,
        reported_proteins=content.reported_proteins,
        raw_files=content.raw_files,
        conditions=content.conditions,
        settings=content.settings,
    )

    hashes = [
        write_result_performance_bundle(
            _ion_intermediate(),
            _proteobot_datapoint(),
            tmp_path / str(index) / "result_performance.csv",
            content=item,
        ).proteobot_json.stem
        for index, item in enumerate((content, reordered, changed_values))
    ]
    assert hashes[0] == hashes[1]
    assert hashes[0] != hashes[2]


def test_content_hash_changes_with_scoring_inputs(tmp_path: Path) -> None:
    content = _submission_content()
    changed_proteins = replace(content, reported_proteins=pd.Series(["P2_HUMAN", "P1_HUMAN"]))
    changed_conditions = replace(content, conditions=("A", "A"))
    changed_settings = replace(content, settings={"species_mapper": {"_HUMAN": "HUMAN"}})
    inputs = (content, changed_proteins, changed_conditions, changed_settings)
    hashes = {
        write_result_performance_bundle(
            _ion_intermediate(),
            _proteobot_datapoint(),
            tmp_path / str(index) / "result_performance.csv",
            content=item,
        ).proteobot_json.stem
        for index, item in enumerate(inputs)
    }

    assert len(hashes) == len(inputs)


def test_content_hash_normalizes_float_precision_nan_and_signed_zero(tmp_path: Path) -> None:
    content = _submission_content()
    float32 = replace(content, matrix=content.matrix.astype(np.float32))
    positive_zero = replace(content, matrix=np.array([[1.0, np.nan], [2.0, 0.0]]))
    hashes = [
        write_result_performance_bundle(
            _ion_intermediate(),
            _proteobot_datapoint(),
            tmp_path / str(index) / "result_performance.csv",
            content=item,
        ).proteobot_json.stem
        for index, item in enumerate((content, float32, positive_zero))
    ]

    assert len(set(hashes)) == 1


def test_bundle_rejects_caller_owned_hash(tmp_path: Path) -> None:
    target = tmp_path / "result_performance.csv"
    datapoint = {**_proteobot_datapoint(), "intermediate_hash": "a" * 40}

    with pytest.raises(ValueError, match="identity fields belong"):
        write_result_performance_bundle(
            _ion_intermediate(), datapoint, target, content=_submission_content()
        )

    assert not target.exists()


def test_write_result_performance_bundle_refuses_either_existing_output(tmp_path: Path) -> None:
    target = tmp_path / "result_performance.csv"
    content = _submission_content()
    seed = write_result_performance_bundle(
        _ion_intermediate(),
        _proteobot_datapoint(),
        tmp_path / "seed" / "result_performance.csv",
        content=content,
    )
    json_target = tmp_path / seed.proteobot_json.name
    json_target.write_text("owned by user\n", encoding="utf-8")

    with pytest.raises(ValueError, match="already exists"):
        write_result_performance_bundle(
            _ion_intermediate(), _proteobot_datapoint(), target, content=content
        )

    assert not target.exists()
    assert json_target.read_text(encoding="utf-8") == "owned by user\n"


def test_write_result_performance_bundle_rolls_back_partial_publication(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    target = tmp_path / "result_performance.csv"
    real_link = os.link
    calls = 0

    def fail_second_link(source: Path, destination: Path) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise FileExistsError(errno.EEXIST, "simulated race", destination)
        real_link(source, destination)

    content = _submission_content()
    seed = write_result_performance_bundle(
        _ion_intermediate(),
        _proteobot_datapoint(),
        tmp_path / "seed" / "result_performance.csv",
        content=content,
    )
    monkeypatch.setattr(os, "link", fail_second_link)
    with pytest.raises(ValueError, match="already exists"):
        write_result_performance_bundle(
            _ion_intermediate(), _proteobot_datapoint(), target, content=content
        )

    assert not target.exists()
    assert not target.with_name(seed.proteobot_json.name).exists()
    assert list(tmp_path.glob(".*.tmp")) == []

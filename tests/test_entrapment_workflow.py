"""Entrapment scoring from vendor files through APB2, apb-catalog and the CLI."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import polars as pl
import pytest
from apb2.api import ParseRuleCompiler, read_parsed_levels

from apb_proteobench.api import EntrapmentAnalyzer
from apb_proteobench.calculation.entrapment import read_pairs
from apb_proteobench.cli import app
from apb_proteobench.configuration.entrapment import load_packaged_entrapment_module
from apb_proteobench.entrapment import CATALOGUE, DIAGNOSTICS_SLOT

COLUMNS = (
    "Run",
    "Modified.Sequence",
    "Stripped.Sequence",
    "Precursor.Charge",
    "Precursor.Id",
    "Protein.Group",
    "Protein.Ids",
    "Protein.Names",
    "Genes",
    "Precursor.Normalised",
    "Q.Value",
    "Lib.Q.Value",
    "Global.Q.Value",
)
# Run, modified sequence, run q-value, library q-value, global q-value. Pair 0: the
# entrapment loses on every kind; pair 1: tied on run and library, better on global;
# pair 2: worse on run, better on library and global. GGMC's carbamidomethyl sits on its
# last residue.
ROWS = (
    ("R1", "AAAK", 0.001, 0.001, 0.001),
    ("R2", "AAAK", 0.003, 0.001, 0.001),
    ("R1", "CCCK", 0.002, 0.002, 0.002),
    ("R1", "DDDK", 0.003, 0.003, 0.003),
    ("R1", "EEEK", 0.003, 0.003, 0.0005),
    ("R1", "GGMC(UniMod:4)", 0.004, 0.004, 0.004),
    ("R1", "HHHK", 0.005, 0.0005, 0.0005),
)
PAIRS = """sequence\tdecoy\tproteins\tpeptide_type\tpeptide_pair_index
AAAK\tNo\tsp|P1|ONE_HUMAN\ttarget\t0
CCCK\tNo\tsp|P1_p_target|ONE_HUMAN_p_target\tp_target\t0
DDDK\tNo\tsp|P2|TWO_HUMAN\ttarget\t1
EEEK\tNo\tsp|P2_p_target|TWO_HUMAN_p_target\tp_target\t1
GGMC\tNo\tsp|P3|THREE_HUMAN\ttarget\t2
GGMC[Carbamidomethyl]\tNo\tsp|P3|THREE_HUMAN\ttarget\t2
HHHK\tNo\tsp|P3_p_target|THREE_HUMAN_p_target\tp_target\t2
"""


def _write_inputs(folder: Path) -> tuple[Path, Path, Path, Path]:
    lines = ["\t".join(COLUMNS)]
    for run, modified, q_value, library, experiment in ROWS:
        stripped = modified.split("(", 1)[0]
        lines.append(
            "\t".join(
                (
                    run,
                    modified,
                    stripped,
                    "2",
                    f"{modified}2",
                    "P1",
                    "P1",
                    "P1_HUMAN",
                    "GENE1",
                    "1000",
                    str(q_value),
                    str(library),
                    str(experiment),
                )
            )
        )
    data = folder / "report.tsv"
    data.write_text("\n".join(lines) + "\n", encoding="utf-8")
    parameters = folder / "report.log.txt"
    parameters.write_text(
        "DIA-NN 1.8.1 (Data-Independent Acquisition by Neural Networks)\ndiann --unimod4\n",
        encoding="utf-8",
    )
    pairs = folder / "pairs.txt.gz"
    pairs.write_bytes(gzip.compress(PAIRS.encode()))
    fasta = folder / "peptides.fasta"
    fasta.write_text(">sp|P1|ONE_HUMAN All peptides\nMAAAKCCCKDDDKEEEKGGMCHHHK\n", encoding="utf-8")
    return data, parameters, pairs, fasta


def test_every_catalogued_q_value_kind_is_scored(tmp_path: Path) -> None:
    data, parameters, pairs, _fasta = _write_inputs(tmp_path)
    parsed = ParseRuleCompiler(data, parameters, requested_levels=("ion",)).compile().parse()
    analyzer = EntrapmentAnalyzer(
        load_packaged_entrapment_module("entrapment_dia_astral"), pairs=read_pairs(pairs)
    )

    result = analyzer.analyze(parsed)

    assert {kind: scores.paired_FDP for kind, scores in result.scores.items()} == {
        "q_value": 3 / 6,
        "library_q_value": 5 / 6,
        "global_q_value": 7 / 6,
    }
    run = result.scores["q_value"]
    assert (run.nr_id_features, run.lower_bound_FDP, run.combined_FDP) == (6, 0.5, 1.0)
    assert run.reported_fdr_parsed_from_input == pytest.approx(0.005)
    sequences = dict(result.precursors.select("peptide", "sequence").iter_rows())
    assert sequences["GGMC"] == "GGMC[Carbamidomethyl]", "last-residue modification"
    best = dict(result.precursors.select("peptide", "q_value").iter_rows())
    assert best["AAAK"] == pytest.approx(0.001), "best run"

    ion = result.parsed.levels["ion"]
    labels = ion.varm[DIAGNOSTICS_SLOT]
    assert isinstance(labels, pl.DataFrame)
    assert labels.height == ion.var.frame.height
    assert sorted(labels.get_column("label").to_list()) == ["entrapment"] * 3 + ["target"] * 3
    stored = ion.metadata["proteobench"]
    assert isinstance(stored, dict)
    entrapment = stored["entrapment"]
    assert isinstance(entrapment, dict)
    assert sorted(entrapment) == ["global_q_value", "library_q_value", "q_value"]
    catalog = result.parsed.metadata["catalog"]
    assert isinstance(catalog, dict)
    assert CATALOGUE in catalog


def test_a_result_scored_once_refuses_a_second_scoring(tmp_path: Path) -> None:
    data, parameters, pairs, _fasta = _write_inputs(tmp_path)
    parsed = ParseRuleCompiler(data, parameters, requested_levels=("ion",)).compile().parse()
    analyzer = EntrapmentAnalyzer(
        load_packaged_entrapment_module("entrapment_dia_astral"), pairs=read_pairs(pairs)
    )
    scored = analyzer.analyze(parsed).parsed

    with pytest.raises(ValueError, match="already exist"):
        analyzer.analyze(scored)


def test_only_packaged_entrapment_modules_load() -> None:
    with pytest.raises(ValueError, match="unknown packaged entrapment module"):
        load_packaged_entrapment_module("dia_astral")


def test_cli_run_entrapment_writes_a_scored_result(tmp_path: Path) -> None:
    data, parameters, pairs, fasta = _write_inputs(tmp_path)
    target = tmp_path / "scored.h5ad"
    timings = tmp_path / "timings"

    status = app(
        [
            "run",
            "entrapment",
            str(data),
            str(fasta),
            "--params",
            str(parameters),
            "--pairs",
            str(pairs),
            "--output",
            str(target),
            "--timings-dir",
            str(timings),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    stored = read_parsed_levels(target).levels["ion"].metadata["proteobench"]
    assert isinstance(stored, dict)
    entrapment = stored["entrapment"]
    assert isinstance(entrapment, dict)
    library = entrapment["library_q_value"]
    assert isinstance(library, dict)
    assert library["paired_FDP"] == 5 / 6
    benchmark = json.loads((timings / "apb-proteobench.benchmark.timings.json").read_text())
    assert [phase["name"] for phase in benchmark["phases"]] == ["load_module", "analyze", "write"]


def test_cli_run_entrapment_needs_the_pair_file(tmp_path: Path) -> None:
    data, parameters, _pairs, fasta = _write_inputs(tmp_path)

    status = app(
        [
            "run",
            "entrapment",
            str(data),
            str(fasta),
            "--params",
            str(parameters),
            "--output",
            str(tmp_path / "scored.h5ad"),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 1
    assert not (tmp_path / "scored.h5ad").exists()

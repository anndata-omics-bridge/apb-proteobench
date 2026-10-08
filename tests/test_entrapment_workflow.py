"""Entrapment scoring from vendor files through APB2, apb-catalog and the CLI."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest
from apb2.api import ParseRuleCompiler, read_parsed_levels
from apb_fasta.api import FastaAnnotator

from apb_proteobench.api import EntrapmentAnalyzer
from apb_proteobench.cli.app import app
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
# pair 2: worse on run, better on library and global. GGMC's and HHHC's carbamidomethyl sits
# on their last residue.
ROWS = (
    ("R1", "AAAK", 0.001, 0.001, 0.001),
    ("R2", "AAAK", 0.003, 0.001, 0.001),
    ("R1", "CCCK", 0.002, 0.002, 0.002),
    ("R1", "DDDK", 0.003, 0.003, 0.003),
    ("R1", "EEEK", 0.003, 0.003, 0.0005),
    ("R1", "GGMC(UniMod:4)", 0.004, 0.004, 0.004),
    ("R1", "HHHC(UniMod:4)", 0.005, 0.0005, 0.0005),
)
# ProteoBench's entrapment FASTA: each target directly before its entrapment.
PEPTIDES = ("AAAK", "CCCK", "DDDK", "EEEK", "GGMC", "HHHC")


def _write_inputs(folder: Path) -> tuple[Path, Path, Path]:
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
    fasta = folder / "peptides.fasta"
    fasta.write_text(
        "".join(
            f">sp|{p}_{kind}|{p}_{kind}\n{p}\n"
            for p, kind in zip(PEPTIDES, ("target", "p_target") * 3, strict=True)
        ),
        encoding="utf-8",
    )
    return data, parameters, fasta


def _fasta(fasta: Path) -> FastaAnnotator:
    return FastaAnnotator.read((fasta,))


def test_every_catalogued_q_value_kind_is_scored(tmp_path: Path) -> None:
    data, parameters, fasta = _write_inputs(tmp_path)
    parsed = ParseRuleCompiler(data, parameters, requested_levels=("ion",)).compile().parse()
    analyzer = EntrapmentAnalyzer(
        load_packaged_entrapment_module("entrapment_dia_astral"), _fasta(fasta)
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
    record = stored["result"]
    assert isinstance(record, dict)
    entrapment = record["entrapment"]
    assert isinstance(entrapment, dict)
    assert sorted(entrapment) == ["global_q_value", "library_q_value", "q_value"]
    assert stored["details"] == [{"slot": "varm", "name": DIAGNOSTICS_SLOT}]
    summary = stored["summary"]
    assert isinstance(summary, list)
    assert {(entry["name"], entry["layer"]) for entry in summary if isinstance(entry, dict)} == {
        (name, kind)
        for name in ("identified_features", "combined_fdp", "paired_fdp")
        for kind in ("global_q_value", "library_q_value", "q_value")
    }
    catalog = result.parsed.metadata["catalog"]
    assert isinstance(catalog, dict)
    assert CATALOGUE in catalog


def test_a_result_scored_once_refuses_a_second_scoring(tmp_path: Path) -> None:
    data, parameters, fasta = _write_inputs(tmp_path)
    parsed = ParseRuleCompiler(data, parameters, requested_levels=("ion",)).compile().parse()
    analyzer = EntrapmentAnalyzer(
        load_packaged_entrapment_module("entrapment_dia_astral"), _fasta(fasta)
    )
    scored = analyzer.analyze(parsed).parsed

    with pytest.raises(ValueError, match="already exist"):
        analyzer.analyze(scored)


def test_only_packaged_entrapment_modules_load() -> None:
    with pytest.raises(ValueError, match="unknown packaged entrapment module"):
        load_packaged_entrapment_module("dia_astral")


def test_cli_run_entrapment_writes_a_scored_result(tmp_path: Path) -> None:
    data, parameters, fasta = _write_inputs(tmp_path)
    target = tmp_path / "scored.h5ad"
    timings = tmp_path / "timings"
    scores = tmp_path / "scores.json"

    status = app(
        [
            "run",
            "entrapment",
            str(data),
            str(fasta),
            "--params",
            str(parameters),
            "--output",
            str(target),
            "--timings-dir",
            str(timings),
            "--scores",
            str(scores),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    datapoints = json.loads(scores.read_text())
    assert sorted(datapoints) == ["global_q_value", "library_q_value", "q_value"]
    assert datapoints["library_q_value"]["paired_FDP"] == 5 / 6
    assert {"software_name", "precursor_mass_tolerance", "fdp_curve"} <= set(datapoints["q_value"])
    stored = read_parsed_levels(target).levels["ion"].metadata["proteobench"]
    assert isinstance(stored, dict)
    record = stored["result"]
    assert isinstance(record, dict)
    entrapment = record["entrapment"]
    assert isinstance(entrapment, dict)
    library = entrapment["library_q_value"]
    assert isinstance(library, dict)
    assert library["paired_FDP"] == 5 / 6
    benchmark = json.loads((timings / "apb-proteobench.benchmark.timings.json").read_text())
    assert [phase["name"] for phase in benchmark["phases"]] == ["load_module", "analyze", "write"]

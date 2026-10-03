"""Direct vendor conversion through APB2's high-level in-memory API."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from apb2.api import read_parsed_levels

from apb_proteobench.cli import app
from conftest import matrix_values, write_module


def _write_benchmark_diann_input(folder: Path) -> tuple[Path, Path, Path, Path]:
    data = folder / "benchmark-report.tsv"
    columns = (
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
        "Fragment.Quant.Raw",
        "Fragment.Correlations",
        "PG.MaxLFQ",
    )
    peptides = ("PEPTIDE", "YEASTPEP", "ECILIPEP", "CONTPEP", "MIXPEP", "NKPEPT")
    proteins = (
        "P1_HUMAN",
        "P2_YEAST",
        "P3_ECOLI",
        "Cont_P4_HUMAN",
        "P5_HUMAN_YEAST",
        "P6_UNKNOWN",
    )
    runs = ("run_A1", "run_A2", "run_B1", "run_B2")
    values = matrix_values().copy()
    values[values <= 0] = 1.0
    rows = ["\t".join(columns)]
    for run_index, run in enumerate(runs):
        for feature_index, (peptide, protein) in enumerate(zip(peptides, proteins, strict=True)):
            value = values[run_index, feature_index]
            rows.append(
                "\t".join(
                    (
                        run,
                        peptide,
                        peptide,
                        "2",
                        f"{peptide}2",
                        protein,
                        protein,
                        protein,
                        f"GENE{feature_index}",
                        str(value),
                        f"{value};{value}",
                        "0.9;0.8",
                        str(value),
                    )
                )
            )
    data.write_text("\n".join(rows) + "\n", encoding="utf-8")
    parameters = folder / "benchmark-report.log.txt"
    parameters.write_text(
        "DIA-NN 1.8.1 (Data-Independent Acquisition by Neural Networks)\ndiann --unimod4\n",
        encoding="utf-8",
    )
    fasta = folder / "proteins.fasta"
    entries = ("P1|P1_HUMAN", "P2|P2_YEAST", "P3|P3_ECOLI", "Cont_P4|P4_HUMAN", "P5|P5_HUMAN")
    sequences = ("MPEPTIDEK", "YEASTPEPK", "ECILIPEPK", "CONTPEPK", "MIXPEPK")
    records = [
        *zip(entries, sequences, strict=True),
        ("P5Y|P5_YEAST", "MIXPEPK"),
        ("P6|P6_UNKNOWN", "NKPEPTK"),
    ]
    fasta.write_text("".join(f">sp|{entry}\n{sequence}\n" for entry, sequence in records))
    module = folder / "module.toml"
    write_module(module)
    return data, parameters, fasta, module


def test_cli_run_can_convert_only_ion_to_h5ad(tmp_path: Path) -> None:
    data, parameters, fasta, module = _write_benchmark_diann_input(tmp_path)
    target = tmp_path / "stored.h5ad"
    result_performance = tmp_path / "pmultiqc" / "result_performance.csv"
    timings_dir = tmp_path / "timings"

    status = app(
        [
            "run",
            "quant",
            str(data),
            str(fasta),
            "--params",
            str(parameters),
            "--module",
            str(module),
            "--software",
            "diann",
            "--level",
            "ion",
            "--output",
            str(target),
            "--result-performance",
            str(result_performance),
            "--timings-dir",
            str(timings_dir),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    assert list(read_parsed_levels(target).levels) == ["ion"]
    assert result_performance.is_file()
    expected = {
        "apb2.convert.timings.json": ("apb2", "convert", ["compile", "read", "parse"]),
        "apb-fasta.verify-peptides.timings.json": (
            "apb-fasta",
            "verify-peptides",
            ["load_database", "verify_peptides"],
        ),
        "apb-proteobench.benchmark.timings.json": (
            "apb-proteobench",
            "benchmark",
            ["load_module", "analyze", "write", "export"],
        ),
    }
    assert {path.name for path in timings_dir.iterdir()} == set(expected)
    for name, (tool, operation, phases) in expected.items():
        document = json.loads((timings_dir / name).read_text(encoding="utf-8"))
        assert document["format"] == "apb-tool-timings"
        assert document["format_version"] == 1
        assert document["tool"] == tool
        assert document["operation"] == operation
        assert [phase["name"] for phase in document["phases"]] == phases
        assert all(phase["seconds"] >= 0 for phase in document["phases"])
    conversion = json.loads((timings_dir / "apb2.convert.timings.json").read_text())
    assert [level["level"] for level in conversion["levels"]] == ["ion"]


def test_cli_run_refuses_existing_timing_file_before_scoring(tmp_path: Path) -> None:
    data, parameters, fasta, module = _write_benchmark_diann_input(tmp_path)
    timings_dir = tmp_path / "timings"
    timings_dir.mkdir()
    (timings_dir / "apb2.convert.timings.json").write_text("old", encoding="utf-8")
    target = tmp_path / "stored.h5ad"

    status = app(
        [
            "run",
            "quant",
            str(data),
            str(fasta),
            "--params",
            str(parameters),
            "--module",
            str(module),
            "--software",
            "diann",
            "--level",
            "ion",
            "--output",
            str(target),
            "--timings-dir",
            str(timings_dir),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 1
    assert not target.exists()
    assert (timings_dir / "apb2.convert.timings.json").read_text() == "old"


@pytest.mark.parametrize("suffix", [".h5mu", ".parquet", ".duckdb"])
def test_cli_run_starts_at_vendor_files_and_writes_final_result(
    suffix: str,
    tmp_path: Path,
) -> None:
    data, parameters, fasta, module = _write_benchmark_diann_input(tmp_path)
    target = tmp_path / f"stored{suffix}"
    result_performance = tmp_path / "pmultiqc" / "result_performance.csv"

    status = app(
        [
            "run",
            "quant",
            str(data),
            str(fasta),
            "--params",
            str(parameters),
            "--module",
            str(module),
            "--software",
            "diann",
            "--output",
            str(target),
            "--result-performance",
            str(result_performance),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    exported = pd.read_csv(result_performance)
    assert "precursor ion" in exported
    assert len(exported) == 3
    proteobot_paths = list(result_performance.parent.glob("*.json"))
    assert len(proteobot_paths) == 1
    proteobot = json.loads(proteobot_paths[0].read_text(encoding="utf-8"))
    assert proteobot["software_name"] == "DIA-NN"
    assert proteobot["software_version"] == "1.8.1"
    assert proteobot_paths[0].stem == proteobot["intermediate_hash"]
    assert proteobot["fixed_mods"] == "C[Carbamidomethyl]"
    restored = read_parsed_levels(target)
    assert "fasta_validation" in restored.levels["ion"].varm
    assert "proteobench:Precursor_Normalised" in restored.levels["ion"].varm

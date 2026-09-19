"""Direct vendor conversion through APB2's high-level in-memory API."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from apb2.result_facade import read_parsed_levels

from apb_proteobench.api import convert_vendor_result, run_vendor_benchmark
from apb_proteobench.cli import app
from conftest import matrix_values, write_module


def _write_diann_input(folder: Path) -> tuple[Path, Path]:
    data = folder / "report.tsv"
    data.write_text(
        "\t".join(
            (
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
        )
        + "\n"
        + "\t".join(
            (
                "run-A",
                "PEPTC(UniMod:4)IDE",
                "PEPTCIDE",
                "2",
                "PEPTCIDE2",
                "P1",
                "P1",
                "Protein 1",
                "GENE1",
                "100",
                "10;20",
                "0.5;0.7",
                "90",
            )
        )
        + "\n"
        + "\t".join(
            (
                "run-B",
                "PEPTC(UniMod:4)IDE",
                "PEPTCIDE",
                "2",
                "PEPTCIDE2",
                "P1",
                "P1",
                "Protein 1",
                "GENE1",
                "200",
                "30;40",
                "0.8;0.9",
                "180",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    parameters = folder / "report.log.txt"
    parameters.write_text(
        "DIA-NN 1.8.1 (Data-Independent Acquisition by Neural Networks)\ndiann --unimod4\n",
        encoding="utf-8",
    )
    return data, parameters


def _write_alphadia_input(folder: Path) -> tuple[Path, Path, Path]:
    matrix = folder / "input_file.tsv"
    matrix.write_text("mod_seq_charge_hash\trun-A\n1\t100\n", encoding="utf-8")
    precursors = folder / "input_file_secondary.tsv"
    precursors.write_text(
        "\t".join(
            (
                "mod_seq_charge_hash",
                "sequence",
                "charge",
                "mods",
                "mod_sites",
                "genes",
                "proteins",
                "pg",
                "pg_master",
                "decoy",
            )
        )
        + "\n"
        + "1\tPEPTIDE\t2\t\t\tGENE1\tP1\tP1\tP1\t0\n",
        encoding="utf-8",
    )
    parameters = folder / "alphadia.log.txt"
    parameters.write_text("0:00:00.0 PROGRESS: version: 1.12.1\n", encoding="utf-8")
    return matrix, precursors, parameters


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
    fasta.write_text(
        ">sp|ALL|ALL All peptides\nMPEPTIDEKYEASTPEPKECILIPEPKCONTPEPKMIXPEPKNKPEPTK\n",
        encoding="utf-8",
    )
    module = folder / "module.toml"
    write_module(module)
    return data, parameters, fasta, module


@pytest.mark.parametrize("suffix", [".h5ad", ".parquet", ".duckdb"])
def test_convert_vendor_result_writes_one_level_with_compiled_parser(
    suffix: str,
    tmp_path: Path,
) -> None:
    data, parameters = _write_diann_input(tmp_path)
    target = tmp_path / "results" / f"ion{suffix}"

    result = convert_vendor_result(
        data,
        parameters,
        target,
        level="ion",
        software="diann",
    )

    restored = read_parsed_levels(target)
    assert result.software == "diann"
    assert result.software_version == "1.8.1"
    assert list(result.parsed.levels) == ["ion"]
    assert list(restored.levels) == ["ion"]
    assert restored.levels["ion"].obs.frame.height == 2
    assert restored.levels["ion"].var.frame.height == 1
    assert restored.uns == {}
    level_provenance = restored.levels["ion"].uns
    assert "rule_json" in level_provenance
    assert "search_parameters" not in level_provenance
    assert "search_parameters_path" not in level_provenance


def test_convert_vendor_result_reads_alphadia_directory_bundle(tmp_path: Path) -> None:
    matrix, precursors, parameters = _write_alphadia_input(tmp_path)
    target = tmp_path / "alphadia.parquet"

    result = convert_vendor_result(
        tmp_path,
        parameters,
        target,
        level="ion",
        software="alphadia",
    )

    assert matrix.is_file()
    assert precursors.is_file()
    ion = result.parsed.levels["ion"]
    assert ion.obs.frame["run"].to_list() == ["run-A"]
    assert ion.layers["Intensity"].values.row(0)[1:] == (100.0,)


@pytest.mark.parametrize("suffix", [".h5mu", ".parquet", ".duckdb"])
def test_run_vendor_benchmark_writes_one_complete_result(
    suffix: str,
    tmp_path: Path,
) -> None:
    data, parameters, fasta, module = _write_benchmark_diann_input(tmp_path)
    target = tmp_path / f"stored{suffix}"

    result = run_vendor_benchmark(
        data,
        parameters,
        (fasta,),
        module,
        target,
        software="diann",
    )

    restored = read_parsed_levels(target)
    assert result.software == "diann"
    assert result.scored.layers["Precursor_Normalised"].analysis.scores.nr_feature == 3
    assert result.fasta_reports.peptide_levels["ion"].unmatched_feature_count == 0
    assert "fasta_validation" in restored.levels["ion"].varm
    assert "proteobench:Precursor_Normalised" in restored.levels["ion"].varm
    assert restored.levels["ion"].obs.frame.get_column("sample_name").to_list() == [
        "A1",
        "A2",
        "B1",
        "B2",
    ]


def test_run_vendor_benchmark_can_convert_only_ion_to_h5ad(tmp_path: Path) -> None:
    data, parameters, fasta, module = _write_benchmark_diann_input(tmp_path)
    target = tmp_path / "stored.h5ad"

    result = run_vendor_benchmark(
        data,
        parameters,
        (fasta,),
        module,
        target,
        level="ion",
        software="diann",
    )

    assert list(result.parsed.levels) == ["ion"]
    assert list(read_parsed_levels(target).levels) == ["ion"]


def test_cli_run_can_convert_only_ion_to_h5ad(tmp_path: Path) -> None:
    data, parameters, fasta, module = _write_benchmark_diann_input(tmp_path)
    target = tmp_path / "stored.h5ad"
    result_performance = tmp_path / "pmultiqc" / "result_performance.csv"

    status = app(
        [
            "run",
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
            "--x",
            "--output",
            str(target),
            "--result-performance",
            str(result_performance),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    assert list(read_parsed_levels(target).levels) == ["ion"]
    assert result_performance.is_file()


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
            "--x",
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

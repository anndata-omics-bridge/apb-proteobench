"""Direct vendor conversion through the APB2 compiler/parser boundary."""

from __future__ import annotations

import json
from pathlib import Path

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


def test_convert_vendor_result_writes_one_level_with_compiled_parser(tmp_path: Path) -> None:
    data, parameters = _write_diann_input(tmp_path)
    target = tmp_path / "results" / "ion.h5ad"

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
    provenance = restored.levels["ion"].uns
    assert provenance["rule_selection_method"] == "software_version"
    assert provenance["search_parameters_path"] == str(parameters)
    search_parameters = provenance["search_parameters"]
    assert isinstance(search_parameters, str)
    assert json.loads(search_parameters)["software_version"] == "1.8.1"


def test_cli_convert_without_level_writes_every_compatible_level(tmp_path: Path) -> None:
    data, parameters = _write_diann_input(tmp_path)
    output = tmp_path / "results" / "all-levels"

    status = app(
        [
            "convert",
            str(data),
            "--params",
            str(parameters),
            "--software",
            "diann",
            "--output",
            str(output),
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    restored = read_parsed_levels(output.with_suffix(".h5mu"))
    assert list(restored.levels) == ["ion", "protein", "fragment"]


def test_cli_convert_requires_parameter_file(tmp_path: Path) -> None:
    data, _parameters = _write_diann_input(tmp_path)

    status = app(
        ["convert", str(data), "ion"],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 1


def test_run_vendor_benchmark_writes_one_complete_h5mu(tmp_path: Path) -> None:
    data, parameters, fasta, module = _write_benchmark_diann_input(tmp_path)
    target = tmp_path / "stored.h5mu"

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
    assert result.scored.analysis.scores.nr_feature == 3
    assert result.fasta_reports.peptide_levels["ion"].unmatched_feature_count == 0
    assert "fasta_validation" in restored.levels["ion"].varm
    assert "proteobench" in restored.levels["ion"].varm
    assert restored.levels["ion"].obs.frame.get_column("sample_name").to_list() == [
        "A1",
        "A2",
        "B1",
        "B2",
    ]


def test_cli_run_starts_at_vendor_files_and_writes_final_h5mu(tmp_path: Path) -> None:
    data, parameters, fasta, module = _write_benchmark_diann_input(tmp_path)
    target = tmp_path / "stored.h5mu"

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
        ],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    restored = read_parsed_levels(target)
    assert "fasta_validation" in restored.levels["ion"].varm
    assert "proteobench" in restored.levels["ion"].varm

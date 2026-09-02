"""Direct vendor conversion through the APB2 compiler/parser boundary."""

from __future__ import annotations

import json
from pathlib import Path

from apb2.result_facade import read_parsed_levels

from apb_proteobench.api import convert_vendor_result
from apb_proteobench.cli import app


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

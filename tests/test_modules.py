"""Packaged ProteoBench module catalogue."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from importlib.resources import as_file, files
from pathlib import Path

import pytest

from apb_proteobench.configuration.entrapment import load_packaged_entrapment_module
from apb_proteobench.configuration.load import (
    PACKAGED_MODULE_NAMES,
    SUPPORTED_MODULE_NAMES,
    load_packaged_module,
)

EXPECTED_SUPPORTED_MODULE_HASHES = {
    "dda_astral": "c4f5ceda352d4dd24e5105a7ee8e4e7f378335d25937cd3527cb8d5cb7e95251",
    "dda_peptidoform": "f9be47e59ba01fa007ec8f8d0a8fac4bfd5094e00b4eb2f505d91799e6027e1d",
    "dda_qexactive": "d0cb8e3a1dc6fcc88af5e77b4638b4827d8f6e59043972b12566e01f8ac88d74",
    "dia_aif": "34f8eaafa688118285398d0ae26517dc5fdf908eb0707c48976069561d398a48",
    "dia_astral": "072f3e31d24bdacfa33c5936c28dd192d137a190c010ace91f27a85d39691a83",
    "dia_diapasef": "0bb73c0722c34c852b4ac956d39b1f738c4be7d2af9b740e861638c4ddb8e063",
    "dia_plasma": "cd8452bf80df674aea07f3250442ed72c8ed33847193c619cf30df106313f462",
    "dia_singlecell": "a64b297599be4503f68ca04b413aa49bd75b84570e1256a51ddbed83d8abf75d",
    "dia_zenotof": "b1514ad1e5d515f6ec620016abb9f1f836e218cdef711bd396b6ac50ad0e1454",
}

EXPECTED_SDRF_HASHES = {
    "dda_astral": "d63efd1e9f0894ce025128c5e4b2e2e4f64ed626968ecf0da0d37db4ee223494",
    "dda_peptidoform": "ac1addb6e6ea42e58f89eb63d4d5e321830c8ca7e15c7a8228e17dbb70849cb5",
    "dda_qexactive": "ac1addb6e6ea42e58f89eb63d4d5e321830c8ca7e15c7a8228e17dbb70849cb5",
    "dia_aif": "fdf5165d3e3ec79442a915193c24931920b9809262a12c466cdb198250783a76",
    "dia_astral": "9979dd3ebbdc04f6820a5b6866a139fb9aa61da95c402820f0f63dd8606ccb0e",
    "dia_diapasef": "b5b9f831b143539cc42dbd13f6c78a6096e069f3cfdf6071d6a99f24c3321171",
    "dia_plasma": "43c9142f635086db339357fa395a46c031cd57a350df0e7a681a9d31f48cf027",
    "dia_singlecell": "36ae5f085ef283eef9798760432bcfdcc12b52a5c3b3bafd0a475ff3d19431e3",
    "dia_zenotof": "08a2a0a964a998b7347a3ec38dd80789437ac24d612ded4f39d15722839ac7e9",
}

EXPECTED_PACKAGED_MODULE_HASHES = {
    **EXPECTED_SUPPORTED_MODULE_HASHES,
    "denovo_dda_hcd": "d53cd86c02228c57cc35aba996e44c4739ce49e4de8fe57bf556a266d18460b7",
    "entrapment_dia_astral": ("c6740cd8f2954b9b137bc870437331c53eb92c09a0e53650a30884cd79692db1"),
}


def test_every_supported_module_is_packaged_validated_and_pinned() -> None:
    assert tuple(EXPECTED_SUPPORTED_MODULE_HASHES) == SUPPORTED_MODULE_NAMES

    loaded = {name: load_packaged_module(name) for name in SUPPORTED_MODULE_NAMES}

    assert {name: module.source.sha256 for name, module in loaded.items()} == (
        EXPECTED_SUPPORTED_MODULE_HASHES
    )
    assert {name: module.source.sdrf.sha256 for name, module in loaded.items()} == (
        EXPECTED_SDRF_HASHES
    )
    assert {module.settings.general.level for module in loaded.values()} == {
        "ion",
        "peptidoform",
    }
    assert all(
        len(module.settings.samples) == module.settings.general.max_nr_observed
        for module in loaded.values()
    )
    assert {module.source.format for module in loaded.values()} == {"proteobench-module-toml-sdrf"}


def test_every_upstream_module_document_is_packaged_and_pinned() -> None:
    assert tuple(EXPECTED_PACKAGED_MODULE_HASHES) == PACKAGED_MODULE_NAMES

    module_resources = files("apb_proteobench.data.modules")
    hashes = {
        name: hashlib.sha256(module_resources.joinpath(f"{name}.toml").read_bytes()).hexdigest()
        for name in PACKAGED_MODULE_NAMES
    }

    assert hashes == EXPECTED_PACKAGED_MODULE_HASHES


def test_packaged_ratios_preserve_plot_colors() -> None:
    ratios = load_packaged_module("dda_qexactive").settings.species_expected_ratio

    assert ratios["YEAST"].color == "#88ccef"


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("dda_qexactive", {"YEAST": 2.0, "ECOLI": 0.25, "HUMAN": 1.0}),
        ("dia_singlecell", {"YEAST": 0.2, "HUMAN": 1.2}),
        ("dia_plasma", {"YEAST": 1 / 3, "ECOLI": 2.0, "HUMAN": 1.0}),
    ],
)
def test_packaged_ratios_are_quotients_of_sdrf_species_quantities(
    name: str,
    expected: dict[str, float],
) -> None:
    settings = load_packaged_module(name).settings

    assert {
        species: ratio.a_vs_b for species, ratio in settings.species_expected_ratio.items()
    } == expected
    assert list(settings.species_mapper) == [f"_{species}" for species in expected]


def test_peptidoform_module_names_the_qexactive_raw_files() -> None:
    peptidoform = load_packaged_module("dda_peptidoform").settings.samples
    qexactive = load_packaged_module("dda_qexactive").settings.samples

    assert [sample.raw_file for sample in peptidoform] == [sample.raw_file for sample in qexactive]
    assert peptidoform[0].raw_file_aliases == ["abundance_A_1", "A_1"]


@pytest.mark.parametrize("name", [*EXPECTED_SDRF_HASHES, "entrapment_dia_astral"])
def test_packaged_sdrf_passes_the_sdrf_proteomics_validator(name: str) -> None:
    resource = files("apb_proteobench.data.modules").joinpath(f"{name}.sdrf.tsv")
    templates = ["-t", "ms-proteomics"]
    if "Data-independent acquisition" in resource.read_text(encoding="utf-8"):
        templates += ["-t", "dia-acquisition"]
    validator = Path(sys.executable).with_name("parse_sdrf")

    with as_file(resource) as path:
        completed = subprocess.run(
            [str(validator), "validate-sdrf", "-s", str(path), *templates, "--use_ols_cache_only"],
            capture_output=True,
            check=False,
            text=True,
        )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Everything seems to be fine" in completed.stdout


@pytest.mark.parametrize(
    "name",
    ["denovo_dda_hcd", "entrapment_dia_astral"],
)
def test_packaged_but_unsupported_module_is_explicit(name: str) -> None:
    with pytest.raises(ValueError, match="is not supported by the quantitative scorer"):
        load_packaged_module(name)


def test_unknown_packaged_module_reports_available_names() -> None:
    with pytest.raises(ValueError, match="unknown packaged ProteoBench module 'missing'"):
        load_packaged_module("missing")


def test_entrapment_module_records_its_sdrf() -> None:
    module = load_packaged_entrapment_module("entrapment_dia_astral")

    assert module.sdrf == "entrapment_dia_astral.sdrf.tsv"
    assert module.sdrf_sha256 == "65b1edb018b1fc2a87ab3dfa0bb2e2cac7aa57a20289f014a3418a0a9cc80104"

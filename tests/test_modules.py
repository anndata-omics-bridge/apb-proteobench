"""Packaged ProteoBench module catalogue."""

from __future__ import annotations

import hashlib
from importlib.resources import files

import pytest

from apb_proteobench.configuration.load import (
    available_modules,
    load_packaged_module,
    packaged_module_names,
)

EXPECTED_SUPPORTED_MODULE_HASHES = {
    "dda_astral": "343b757e1c314bdc9dcbdbc3576025d4789fa06779475cfec723e9bae94d1e7f",
    "dda_peptidoform": "62871dfbd14e88cd7f2e98ca74924c65300a0d8198ba13d54d0b9f9ab5c14645",
    "dda_qexactive": "2589e844051bdb43e14a79df64ffdbab7fe1557e75a4298952749b394dd17b05",
    "dia_aif": "5e0b804bc6d979291814338b3274e113b03d2377c88da13161631c3b816c9e8d",
    "dia_astral": "2a78fe4ded14ea7bfe3b386cd900c87a2d4de37d57f08eb6f38c92d52eb49644",
    "dia_diapasef": "2b88d4698b293b77e9a686176442d200265796a7a802372895eaf963149ad221",
    "dia_singlecell": "f695e4f15df838c2b4760079f6afeb4baf8e584c4ace972a81bfb7f2161ae580",
    "dia_zenotof": "635e5610e84ff16241c216be9ade043eb79f03617b60dc8bc9806904b83ed5ed",
}

EXPECTED_PACKAGED_MODULE_HASHES = {
    **EXPECTED_SUPPORTED_MODULE_HASHES,
    "dia_plasma": "11bb9eff64cbbadff19b565bdd1b212b58d86aaeb091f9c9b5ad1493327fdb9d",
    "denovo_dda_hcd": "d53cd86c02228c57cc35aba996e44c4739ce49e4de8fe57bf556a266d18460b7",
    "entrapment_dia_astral": ("dfaa17b3474d91a69eb159ff6a05c59bc89839ab05791f4e8539984d6529f0f2"),
}


def test_every_supported_module_is_packaged_validated_and_pinned() -> None:
    assert available_modules() == tuple(EXPECTED_SUPPORTED_MODULE_HASHES)

    loaded = {name: load_packaged_module(name) for name in available_modules()}

    assert {name: module.source.sha256 for name, module in loaded.items()} == (
        EXPECTED_SUPPORTED_MODULE_HASHES
    )
    assert {module.settings.general.level for module in loaded.values()} == {
        "ion",
        "peptidoform",
    }
    assert all(len(module.settings.samples) == 6 for module in loaded.values())


def test_every_upstream_module_document_is_packaged_and_pinned() -> None:
    assert packaged_module_names() == tuple(EXPECTED_PACKAGED_MODULE_HASHES)

    module_resources = files("apb_proteobench.data.modules")
    hashes = {
        name: hashlib.sha256(module_resources.joinpath(f"{name}.toml").read_bytes()).hexdigest()
        for name in packaged_module_names()
    }

    assert hashes == EXPECTED_PACKAGED_MODULE_HASHES


@pytest.mark.parametrize(
    "name",
    ["dia_plasma", "denovo_dda_hcd", "entrapment_dia_astral"],
)
def test_packaged_but_unsupported_module_is_explicit(name: str) -> None:
    with pytest.raises(ValueError, match="is not supported by the quantitative scorer"):
        load_packaged_module(name)


def test_unknown_packaged_module_reports_available_names() -> None:
    with pytest.raises(ValueError, match="unknown packaged ProteoBench module 'missing'"):
        load_packaged_module("missing")

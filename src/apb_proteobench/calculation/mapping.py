"""ProteoBench-compatible rendering of APB feature identifiers."""

from __future__ import annotations

import re

import pandas as pd
from apb2.api import canonical_modification_names

_UNIMOD_TAG = re.compile(r"\[(UNIMOD:\d+)\]", flags=re.IGNORECASE)
_FINAL_RESIDUE_MODS = re.compile(
    r"(?<=[A-Z])-?(?:\[UNIMOD:\d+\])+(?=(?:/\d+)?$)",
    re.IGNORECASE,
)


def render_proteobench_features(
    features: pd.Series,
    *,
    drop_final_residue_modifications: bool = False,
) -> pd.Series:
    """Render canonical APB ProForma tags with ProteoBench's legacy names.

    APB deliberately retains Unimod accessions in its canonical feature axis,
    whereas ProteoBench's intermediate CSV uses modification names. This
    compatibility rendering is used only for the reconstructed intermediate;
    it never changes the AnnData feature identifiers.
    """
    names = {accession.upper(): name for accession, name in canonical_modification_names().items()}

    def replace(match: re.Match[str]) -> str:
        accession = match.group(1).upper()
        name = names.get(accession)
        return f"[{name}]" if name is not None else match.group(0)

    rendered = features.astype("string")
    if rendered.isna().any():
        raise ValueError("ProteoBench feature identifiers must not be missing")
    if drop_final_residue_modifications:
        # ProteoBench 0.17's ``before_aa = false`` parser does not visit the
        # position after the final residue. Preserve that behavior in the
        # compatibility table without changing APB's canonical feature axis.
        rendered = rendered.str.replace(_FINAL_RESIDUE_MODS, "", regex=True)
    return rendered.str.replace(_UNIMOD_TAG, replace, regex=True)

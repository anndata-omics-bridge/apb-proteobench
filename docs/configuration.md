# Module configuration

The authored configuration is the existing ProteoBench module TOML:

```toml
[species_expected_ratio.YEAST]
A_vs_B = 2.0

[species_expected_ratio.HUMAN]
A_vs_B = 1.0

[species_mapper]
_YEAST = "YEAST"
_HUMAN = "HUMAN"

[general]
min_count_multispec = 1
level = "ion"
default_cutoff_min_feature = 1
max_nr_observed = 6

[[samples]]
raw_file = "run_A1"
raw_file_alias = "run_A1_uncalibrated"
sample_name = "A1"
condition = "A"
```

- `species_expected_ratio` maps each reported species to its expected A/B abundance ratio.
- `species_mapper` maps protein-identifier patterns to those same species names.
- `min_count_multispec` controls exclusion of features assigned to multiple species.
- `level` selects the single APB2 quantification level to annotate and score.
- `default_cutoff_min_feature` selects the aggregate projection shown at the top of the score.
- `max_nr_observed` controls the complete cutoff-indexed score range.
- Each `samples` record supplies a matching `raw_file`, optional alias, displayed `sample_name`,
  and condition. Conditions A and B are required.

HY omits ECOLI from both species sections; HYE includes it. The quantitative scorer supports the
eight HYE/HY modules used by the legacy APB integration:

| Name | Acquisition and instrument | Level |
| --- | --- | --- |
| `dda_astral` | DDA Astral | ion |
| `dda_peptidoform` | DDA | peptidoform |
| `dda_qexactive` | DDA Q Exactive | ion |
| `dia_aif` | DIA AIF | ion |
| `dia_astral` | DIA Astral | ion |
| `dia_diapasef` | DIA diaPASEF | ion |
| `dia_singlecell` | DIA low-input/single-cell | ion |
| `dia_zenotof` | DIA ZenoTOF | ion |

List and load them without locating package files:

```python
from apb_proteobench.configuration.load import available_modules, load_packaged_module

print(available_modules())
module = load_packaged_module("dia_singlecell")
```

Three newer ProteoBench module documents are also packaged so this repository owns the complete
upstream catalogue:

| Name | Module | Current status |
| --- | --- | --- |
| `dia_plasma` | DIA plasma/PYE | scoring not implemented |
| `denovo_dda_hcd` | de novo DDA HCD | different schema; not implemented |
| `entrapment_dia_astral` | DIA Astral entrapment | different schema; not implemented |

Use `packaged_module_names()` to inventory all 11 documents. `available_modules()` deliberately
returns only the eight modules accepted by `load_packaged_module()` and the current quantitative
scorer:

```python
from apb_proteobench.configuration.load import packaged_module_names

print(packaged_module_names())
```

The TOMLs are copied from ProteoBench and pinned by checksum. Their source revision, original paths,
and support status are recorded in the packaged `MODULES_NOTICE.md`.

Unknown upstream module sections are tolerated because the source document is shared with
ProteoBench, but only the validated scoring and sample subset is embedded into APB metadata.
The embedded sample design is normalized as a column mapping (`raw_file`, `sample_name`, and
`condition` arrays), rather than TOML's list of tables, so the same JSON-compatible value also
round-trips through AnnData's HDF5 representation. Authored aliases remain traceable to the
checksum-identified source module but are unnecessary after the match has been applied.

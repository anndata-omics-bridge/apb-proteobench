# ProteoBench module settings

The TOML files in `modules/` are copied from the ProteoBench project at revision
`b69dbaa89be332d644e37e3ff225994aab5947df` and are distributed under the Apache License 2.0.

The package contains all 11 module documents present at that revision. Legacy APB downloaded the
eight HYE/HY documents marked as supported. The three newer documents are retained here as
authoritative inputs for planned support; their schemas or calculations are not yet implemented by
the quantitative scorer.

| Packaged name | ProteoBench source path below `proteobench/io/parsing/io_parse_settings/` | Status |
| --- | --- | --- |
| `dda_astral` | `Quant/lfq/DDA/ion/Astral/module_settings.toml` | Supported |
| `dda_peptidoform` | `Quant/lfq/DDA/peptidoform/module_settings.toml` | Supported |
| `dda_qexactive` | `Quant/lfq/DDA/ion/QExactive/module_settings.toml` | Supported |
| `dia_aif` | `Quant/lfq/DIA/ion/AIF/module_settings.toml` | Supported |
| `dia_astral` | `Quant/lfq/DIA/ion/Astral/module_settings.toml` | Supported |
| `dia_diapasef` | `Quant/lfq/DIA/ion/diaPASEF/module_settings.toml` | Supported |
| `dia_singlecell` | `Quant/lfq/DIA/ion/lowinput/module_settings.toml` | Supported |
| `dia_zenotof` | `Quant/lfq/DIA/ion/ZenoTOF/module_settings.toml` | Supported |
| `dia_plasma` | `Quant/lfq/DIA/ion/plasma/module_settings.toml` | Packaged; scoring planned |
| `denovo_dda_hcd` | `denovo/DDA/HCD/module_settings.toml` | Packaged; schema unsupported |
| `entrapment_dia_astral` | `entrapment/DIA/ion/Astral/module_settings.toml` | Packaged; schema unsupported |

Upstream project: <https://github.com/Proteobench/ProteoBench>

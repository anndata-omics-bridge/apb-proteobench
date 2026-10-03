# ProteoBench module settings

The TOML files in `modules/` are copied from the ProteoBench project at revision
`b69dbaa89be332d644e37e3ff225994aab5947df` and are distributed under the Apache License 2.0.

The package contains all 11 module documents present at that revision. Legacy APB downloaded the
eight HYE/HY documents; together with the plasma document they are marked as supported. The two
other documents are retained here as authoritative inputs for planned support; their schemas are
not yet implemented by the quantitative scorer.

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
| `dia_plasma` | `Quant/lfq/DIA/ion/plasma/module_settings.toml` | Supported |
| `denovo_dda_hcd` | `denovo/DDA/HCD/module_settings.toml` | Packaged; schema unsupported |
| `entrapment_dia_astral` | `entrapment/DIA/ion/Astral/module_settings.toml` | Supported by the entrapment scorer |

On 2026-10-02 the eight supported documents were reshaped for APB: their `[[samples]]` design and species quantities moved into an adjacent `<name>.sdrf.tsv`, species ratios are computed from those quantities, and the TOML keeps the scoring settings plus `[run_aliases]`. Instrument, acquisition, enzyme and mixture facts in the SDRFs come from the ProteoBench module documentation, PRIDE PXD028735, PXD049412 and PXD070049, and the PSI-MS ontology. `dia_plasma` followed on 2026-10-02: its SDRF records PYE9 of Distler et al., Nat. Commun. 2025 (jPOST JPST003358), with instrument, load, enzyme and mixture taken from that paper's Supplementary Information; its HYE FASTA URL follows the current upstream module settings. `entrapment_dia_astral` gained an SDRF for its three 50 ng HeLa runs on an Orbitrap Astral (ProteoBench module documentation); its TOML otherwise stays the upstream copy. `denovo_dda_hcd` remains an unchanged upstream copy.

Upstream project: <https://github.com/Proteobench/ProteoBench>

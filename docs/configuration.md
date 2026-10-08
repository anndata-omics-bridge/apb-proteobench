# Module configuration

## Authored documents

A module is a TOML file and the SDRF-Proteomics table it names. The SDRF owns the sample design: one row per run, its data file, assay name, condition, and the quantity of every species in the mixture. The TOML holds only the scoring settings an SDRF cannot express.

```toml
sdrf = "module.sdrf.tsv"

[species.YEAST]
organism = "saccharomyces cerevisiae"
suffix = "_YEAST"
color = "#88ccef"

[species.HUMAN]
organism = "homo sapiens"
suffix = "_HUMAN"

[general]
min_count_multispec = 1
level = "ion"
default_cutoff_min_feature = 1
max_nr_observed = 6

[run_aliases]
"run_A1.raw" = ["run_A1_uncalibrated"]
```

The SDRF path is relative to the TOML's directory. APB2 reads it through `apb2.api.SdrfSource.read`; the columns ProteoBench uses are:

| SDRF column | Becomes |
| --- | --- |
| `comment[data file]` | run identifier `raw_file`, without directories or extension |
| `assay name` | `sample_name` |
| the single `factor value[…]` column | `condition`; A and B are required |
| `characteristics[spiked compound]`, one column per species | `CT=mixture;SP=<organism>;QY=<amount> <unit>` |

## Fields used by scoring

- Each species' expected A/B ratio is the quotient of its condition A and condition B quantities. Every scored species needs one quantity per condition, in one unit, identical across that condition's rows.
- `[species.<NAME>]` links an SDRF organism to the FASTA entry-name `suffix` that marks its proteins, and to an optional plot `color`. Every spiked organism needs an entry.
- Species and contaminants come from apb-fasta's `varm["fasta_validation"]`, so scoring requires peptide verification. A feature belongs to a species when a FASTA protein containing its peptide has that organism (`suffix` without `_`, e.g. `HUMAN`). It is a contaminant when any such protein is, or when apb2 marked it `apb_Contaminant`; it is dropped as a decoy when apb2 marked it `apb_Decoy`.
- `[run_aliases]` lists further run names that vendor tables report, keyed by SDRF data file. An observation may match a run's identifier, an alias, or its sample name; successful matching always records the module's canonical metadata.
- `min_count_multispec` controls exclusion of features assigned to multiple species.
- `level` selects the single APB2 quantification level to annotate and score.
- `default_cutoff_min_feature` selects the aggregate projection shown at the top of the score.
- `max_nr_observed` controls the complete cutoff-indexed score range.

Unknown TOML sections, such as ProteoBench's `[reference_database]` and `[validation]`, are accepted. Only the resolved sample design and scoring fields are embedded in APB metadata. The packaged SDRFs pass the official `sdrf-pipelines` structural validation; APB itself does not validate SDRF structure or ontology terms.

## Supported quantitative modules

HY omits ECOLI from the SDRF and the `[species]` table; HYE includes it. The quantitative scorer supports the
eight HYE/HY modules used by the legacy APB integration and the plasma module:

| Name | Acquisition and instrument | Level |
| --- | --- | --- |
| `dda_astral` | DDA Astral | ion |
| `dda_peptidoform` | DDA | peptidoform |
| `dda_qexactive` | DDA Q Exactive | ion |
| `dia_aif` | DIA AIF | ion |
| `dia_astral` | DIA Astral | ion |
| `dia_diapasef` | DIA diaPASEF | ion |
| `dia_plasma` | DIA timsTOF, human plasma background (PYE) | ion |
| `dia_singlecell` | DIA low-input/single-cell | ion |
| `dia_zenotof` | DIA ZenoTOF | ion |

### Load from Python

Load them by name without locating package files:

```python
from apb_proteobench.api import load_packaged_module

module = load_packaged_module("dia_singlecell")
```

`load_packaged_module()` accepts only these nine modules, which are validated for the current quantitative scorer; any other name raises and lists them.

## Packaged for planned support

Two newer ProteoBench module documents are also packaged so this repository owns the complete
upstream catalogue:

| Name | Module | Current status |
| --- | --- | --- |
| `denovo_dda_hcd` | de novo DDA HCD | different schema; not implemented |
| `entrapment_dia_astral` | DIA Astral entrapment | scored by `run entrapment`; TOML plus SDRF of the three HeLa runs |

`load_packaged_module()` refuses both. `load_packaged_entrapment_module("entrapment_dia_astral")` loads the entrapment module for `EntrapmentAnalyzer`; see the [Python API](api.md#score-entrapment).

The nine supported modules are derived from ProteoBench's module TOMLs, with their sample design moved into packaged SDRFs; the two others are unchanged copies. All files are pinned by checksum. Their source revision, original paths,
and support status are recorded in
[`MODULES_NOTICE.md`](https://github.com/anndata-omics-bridge/apb-proteobench/blob/main/src/apb_proteobench/data/MODULES_NOTICE.md).

## Persisted form

The embedded sample design is normalized as a column mapping (`raw_file`, `sample_name`, and
`condition` arrays), rather than TOML's list of tables, so the same JSON-compatible value
round-trips through AnnData's HDF5 representation. Authored aliases participate in annotation and
remain traceable to the checksum-identified source module; `source.sdrf` records the SDRF's name and checksum beside the TOML's.

Continue to [Result layout](results.md) for the physical and storage-neutral metadata locations.

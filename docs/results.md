# Result layout

`apb-proteobench` preserves the APB2 result model and adds only namespaced annotation and benchmark
data. The same logical values round-trip through h5ad, h5mu, APB2 Parquet, and DuckDB.

## After annotation

The level selected by `general.level` gains three observation columns:

| Column | Meaning |
| --- | --- |
| `raw_file` | module raw-file identifier used for the match |
| `sample_name` | displayed ProteoBench sample name |
| `condition` | experimental condition used by scoring |

The normalized module configuration and a single source descriptor, including the authored TOML's SHA-256, live at root `apb.proteobench.provenance.annotation`, beside the record's `schema_version` 4. Sample matching runs through APB2's prolfquapp convention, so each level's `apb.prolfquapp` record retains matching counts, corrections, added-column names and their summary. In H5MU the root is MuData; a standalone H5AD keeps the same root part in `uns["apb"]` and the level part in `uns[<level>]["apb"]`.

## After scoring

Scoring writes the `proteobench` record, schema version 4, with common provenance on the root and layer-specific results on the configured level:

| Location | Contents |
| --- | --- |
| `varm["proteobench:<layer-name>"]` | feature-aligned mixed-species diagnostics for one selected layer |
| Root `uns["apb"]["proteobench"]["provenance"]["scoring"]` | versions, method identities, selection mode |
| Level `["proteobench"]["result"]["scoring"]["Intensity"]` | retained quantity name, roles and scores |
| Level `["proteobench"]["summary"]` | per layer: features scored, median and mean absolute epsilon; per q-value kind after entrapment: identified precursors, combined and paired FDP |
| Level `["proteobench"]["details"]` | the `varm` slots holding each layer's diagnostics |

The feature diagnostics have exactly one row per variable in the selected level. Scores include the
complete cutoff-indexed result plus the configured default-cutoff projection.

Each score document follows ProteoBench 0.18 for every quantitative module, including the plasma metrics: per cutoff, spike-in error (`*_abs_log2_fc_error_spike_ins`, plus `_global` and `_eq_species`), counts (`nr_quantified_spike_ins`, `nr_quantified_<SPECIES>`), HUMAN dynamic range (`dynamic_range_human_plasma_A`, `_B`, `_mean`: P90 minus P10 of log10 mean intensity) and HUMAN error (`*_abs_epsilon_human_plasma`). Spike-ins are all species other than HUMAN; outside `dia_plasma`, `human_plasma` names the module's human sample. The default cutoff projects `median_abs_log2_fc_error_spike_ins`, `nr_quantified_spike_ins`, `dynamic_range_human_plasma` and `median_abs_epsilon_human_plasma` to the top level. As in ProteoBench, an empty row subset scores 0.0. Where APB's inputs differ from ProteoBench's, see [Differences from ProteoBench](compatibility.md).

Annotation remains intact when scoring is added. Scoring refuses existing scores, not an annotation-only namespace. Selected quantities are derived from scoring keys: there is no `layers` wrapper or `X` alias. The primary quantity retains its logical name, displayed as `Intensity · X`; its matrix is stored only in `X`.

## pMultiQC compatibility export

`--result-performance reports/result_performance.csv` writes an external compatibility bundle for the selected ion-level layer: the existing `analysis.diagnostics.legacy` table for pMultiQC and a sibling `<intermediate_hash>.json` ProteoBot datapoint matching the `Proteobench/Results_quant_ion_DDA` filename and field conventions. The CSV preserves column order and missing values without adding an index; the writer computes a versioned SHA-256 hash of the selected layer and scoring inputs only for this export. The JSON combines that hash and the completed scores with APB search-parameter metadata. New hashes do not match historical SHA-1 uploads. Neither file changes the APB2 storage schema or adds a pMultiQC dependency. The writer refuses non-canonical CSV filenames or either existing target, stages both files before publication, creates each final path atomically, and rolls back a partial publication on failure.

## Inspect with Python

Read any supported result through APB2's storage-neutral facade:

```python
from pathlib import Path

from apb2.api import read_parsed_levels

parsed = read_parsed_levels(Path("results/scored.h5mu"))
level = parsed.levels["ion"]

diagnostics = level.varm["proteobench:Intensity"]
score_record = level.metadata["proteobench"]
intensity_record = score_record["result"]["scoring"]["Intensity"]

print(diagnostics.head())
print(parsed.metadata["proteobench"]["provenance"]["scoring"])
print(intensity_record["scores"])
```

The in-memory model keeps root provenance in `ParsedLevels.metadata["proteobench"]["provenance"]` and level results in `ParsedLevel.metadata["proteobench"]`. Metadata scoring keys retain reversible percent encoding for backend-unsafe characters, while each record's `layer_name` retains the exact logical name. APB2 maps diagnostic slot names to each physical backend without requiring AnnData or MuData in the scientific calculation. This is a coordinated pre-1.0 breaking layout change; older persisted formats are rejected without compatibility aliases or automatic migration.

## Supported result formats

| Format | Path | Levels | ProteoBench behavior |
| --- | --- | --- | --- |
| AnnData | `.h5ad` | exactly one | annotation, diagnostics, and scores on the single level |
| MuData | `.h5mu` | one or more | only the module-selected modality is annotated and scored |
| APB2 Parquet | `.parquet` directory | one or more | logical tables and metadata retained exactly |
| DuckDB | `.duckdb` file | one or more | logical tables and metadata retained exactly |

These readers consume APB2-authored results; they are not general importers for arbitrary AnnData,
MuData, Parquet, or DuckDB files. Use APB2 vendor conversion to create the initial result.

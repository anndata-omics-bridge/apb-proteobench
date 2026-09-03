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

The APB annotation namespace retains match provenance, the normalized module configuration, and a
SHA-256 identity of the authored TOML. In h5ad and h5mu this extension metadata is represented below
`uns["apb"]["annotation"]["proteobench"]`.

## After scoring

Scoring adds two values to the configured level:

| Location | Contents |
| --- | --- |
| `varm["proteobench"]` | feature-aligned mixed-species diagnostics |
| `uns["apb"]["proteobench"]` | aggregate scores, methods, column roles, compatibility version, source revision, and protein-mapping provenance |

The feature diagnostics have exactly one row per variable in the selected level. Scores include the
complete cutoff-indexed result plus the configured default-cutoff projection.

## Inspect with Python

Read any supported result through APB2's storage-neutral facade:

```python
from pathlib import Path

from apb2.result_facade import read_parsed_levels

parsed = read_parsed_levels(Path("results/scored.h5mu"))
level = parsed.levels["ion"]

diagnostics = level.varm["proteobench"]
score_record = level.metadata["proteobench"]

print(diagnostics.head())
print(score_record["scores"])
```

The in-memory `ParsedLevels` representation uses `ParsedLevels.metadata["annotation"]` for shared
annotation evidence and `ParsedLevel.metadata["proteobench"]` for the level-specific score record.
The APB2 writers map those values to each backend without requiring AnnData or MuData in the
scientific calculation.

## Supported result formats

| Format | Path | Levels | ProteoBench behavior |
| --- | --- | --- | --- |
| AnnData | `.h5ad` | exactly one | annotation, diagnostics, and scores on the single level |
| MuData | `.h5mu` | one or more | only the module-selected modality is annotated and scored |
| APB2 Parquet | `.parquet` directory | one or more | logical tables and metadata retained exactly |
| DuckDB | `.duckdb` file | one or more | logical tables and metadata retained exactly |

These readers consume APB2-authored results; they are not general importers for arbitrary AnnData,
MuData, Parquet, or DuckDB files. Use APB2 vendor conversion to create the initial result.

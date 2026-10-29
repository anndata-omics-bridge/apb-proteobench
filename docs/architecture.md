# Architecture

## Dependency direction

Dependency direction follows ownership:

```text
CLI / presentation
        ↓
in-memory API and APB integration
        ↓
calculation protocols and workflow
        ↓
diagnostics / scoring / configuration

apb-proteobench → apb2.api, apb_fasta.api, apb_catalog.api
APB2, APB FASTA, APB Catalog -/→ apb-proteobench
```

## Package boundaries

APB2 owns storage-neutral `ParsedLevels`, result adapters, and generic relational observation
annotation. This package owns module semantics, complete ProteoBench coverage, species mapping,
diagnostics, scoring, compatibility assets, and presentation.

The public API owns one operation: annotate and score canonical `ParsedLevels` in memory. It contains no paths, APB2 compiler, FASTA loading, persistence, logging, or command behavior. The CLI is the physical composition root: raw-vendor workflows construct APB2's `ParseRuleCompiler`, pass canonical parser output to APB FASTA and `ProteoBenchAnalyzer`, persist the returned `ParsedLevels` through APB2, and use `compiler.parameters` for ProteoBot fields without provenance reconstruction. The integration boundary alone extracts scientific input from `ParsedLevels` or attaches calculation output; presentation receives the completed analysis and explicit command paths and never reopens output files.

## Scientific input

Scientific methods receive `QuantitativeLevelInput`: a pandas observation table, NumPy/SciPy
matrix, feature identifiers, protein assignments, and level name. They never receive AnnData,
MuData, or `ParsedLevels`. The integration boundary extracts those values and attaches the completed
analysis to a new `ParsedLevels`; APB2 persists it when the caller requests storage.

## Calculation methods

`DiagnosticMethod` and `ScoringMethod` are internal extension points in `workflow.py`, not part of `api.py`. The default composition uses `MixedSpeciesDiagnostics` and `ProteoBenchCompatibleScoring` and selects no method by name. Automatic plugin discovery is not part of this release.

## Persisted extension data

ProteoBench owns one tool namespace. Root `provenance` records annotation configuration and common scoring methods; the module-selected level records `result.scoring[quantity_name]`, a `summary` and the `details` naming its `varm` slots, while sample matching evidence sits in the level's `prolfquapp` record. Aligned diagnostics remain in that level's `varm`. APB2's writers map these logical owners to H5AD, H5MU, Parquet or DuckDB without tool-specific branches; calculations remain storage-independent. See [Result layout](results.md).

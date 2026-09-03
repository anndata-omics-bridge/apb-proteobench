# Architecture

## Dependency direction

Dependency direction follows ownership:

```text
CLI / presentation
        ↓
API and APB integration
        ↓
calculation protocols and workflow
        ↓
diagnostics / scoring / configuration

apb-proteobench → APB2 public annotation and result facades
APB2 -/→ apb-proteobench
```

## Package boundaries

APB2 owns storage-neutral `ParsedLevels`, result adapters, and generic relational observation
annotation. This package owns module semantics, complete ProteoBench coverage, species mapping,
diagnostics, scoring, compatibility assets, and presentation.

The file-to-file API owns physical orchestration. The integration boundary is the only layer that
extracts scientific input from `ParsedLevels` or attaches calculation output. The CLI composes that
API and presents results; it does not reopen output files.

## Scientific input

Scientific methods receive `QuantitativeLevelInput`: a pandas observation table, NumPy/SciPy
matrix, feature identifiers, protein assignments, and level name. They never receive AnnData,
MuData, or `ParsedLevels`. The integration boundary extracts those values and later persists the
typed result.

## Replaceable methods

`DiagnosticMethod` and `ScoringMethod` are client-owned protocols. The default composition uses
`MixedSpeciesDiagnostics` and `ProteoBenchCompatibleScoring`; callers can explicitly substitute
another implementation without adding method-name branches. Automatic plugin discovery is not
part of this release.

## Persisted extension data

Annotation evidence is shared result metadata. Diagnostics and scores belong to the module-selected
level. APB2's result writers map those logical values to h5ad, h5mu, Parquet, or DuckDB; the
calculation remains independent of every storage backend. See [Result layout](results.md).

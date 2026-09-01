# Architecture

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

APB2 owns storage-neutral `ParsedLevels`, result adapters, and generic relational observation
annotation. This package owns module semantics, complete ProteoBench coverage, species mapping,
diagnostics, scoring, compatibility assets, and presentation.

Scientific methods receive `QuantitativeLevelInput`: a pandas observation table, NumPy/SciPy
matrix, feature identifiers, protein assignments, and level name. They never receive AnnData,
MuData, or `ParsedLevels`. The integration boundary extracts those values and later persists the
typed result.

`DiagnosticMethod` and `ScoringMethod` are client-owned protocols. The default composition uses
`MixedSpeciesDiagnostics` and `ProteoBenchCompatibleScoring`; callers can explicitly substitute
another implementation without adding method-name branches. Automatic plugin discovery is not
part of this release.

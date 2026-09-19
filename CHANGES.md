# Changes

## Unreleased

- **Breaking:** tool-owned metadata combines annotation and scoring under `proteobench`. Root `provenance.annotation` stores configuration/source once; root `provenance.scoring` stores common versions, methods and selection mode (schema 3). Level `annotation` preserves matching evidence and `scoring[quantity_name]` retains per-layer results without a `layers` wrapper or `X` alias. Annotation-only results remain scoreable; existing scores are still protected. Calculations, CLI options, CSV/ProteoBot JSON exports and pMultiQC behavior are unchanged.

- Adapted conversion and scoring to APB2's canonical shared/level metadata scopes; shared search-parameter provenance is no longer copied into every level.
- **Breaking:** Scoring now includes every declared abundance layer by default. The CLI replaces redundant `--all` with explicit `--x` for ProteoBench-compatible primary/X-only scoring; the Python API accepts `PRIMARY_LAYER` for the same selection. Compatibility export requires `--x` or one `--layer NAME`.
- `run` now accepts `--level` to convert one quantification level before FASTA verification and scoring. A one-level run can persist directly as H5AD; omitting the option retains all compatible levels.
- Module sample annotation now accepts `raw_file`, `raw_file_alias`, or `sample_name` as the observed run identifier while preserving the module's `raw_file` and `sample_name` as canonical metadata.
- Added an explicit ion-level export bundle for the unchanged pMultiQC ProteoBench module and ProteoBot result repository. The `run` and `benchmark` commands accept `--result-performance` for one selected layer and publish `result_performance.csv` without an index plus the matching `<intermediate_hash>.json` ProteoBot datapoint. Both files are staged before publication, each final path is created atomically without overwrite, and a failed publication rolls back files added by the same call.
- **Breaking:** Reduced the CLI to the complete `run` and `benchmark` workflows. Staged shell workflows use the dedicated APB commands.
- **Breaking:** Replaced the path-oriented Python functions and their `ConvertedVendorResult`, `ScoredResult`, and `VendorBenchmarkResult` wrappers with `ProteoBenchAnalyzer`. The analyzer binds one validated module and calculation collaborators, accepts canonical `ParsedLevels`, and returns an in-memory `ProteoBenchAnalysisResult` containing only configuration, resolved selection, layer analyses, and the scored APB artifact. APB2 conversion, FASTA verification, persistence, paths, search parameters, logging, and CLI evidence remain at their owning composition boundaries.
- **Breaking:** ProteoBench analysis now returns ordered `ProteoBenchAnalysisResult.layers` entries and persists schema version 3 with per-layer diagnostics at `varm["proteobench:<layer-name>"]` plus layer-keyed records under `metadata["proteobench"]["scoring"]`; the singular result properties and singular `varm["proteobench"]` layout were removed without compatibility aliases.
- The raw-vendor `run` CLI keeps conversion, FASTA annotation, and scoring as storage-neutral APB2 `ParsedLevels`, then persists once through `write_parsed_levels`. Final H5AD, H5MU, Parquet, and DuckDB targets are supported; ProteoBench no longer imposes an H5MU handoff.
- Calculations always compute in float64. A value-dependent heuristic previously inferred the
  arithmetic precision by testing whether the data survived a float32 round-trip, so a float64
  matrix of round numbers silently lost half its mantissa. Storage precision is not arithmetic
  precision: callers may still pass float32, which widens on entry. Removed `FloatDType`,
  `_is_float32_backed`, `_as_float_array` and the `source_dtype` threading.
  On the single-cell HY fixture the epsilon residual falls from 1.8e-07 to 1.1e-16, one ulp,
  which also removes an x86_64/arm64 platform difference that exceeded the test tolerance.
  The legacy golden intermediate still matches within 7.1e-08 -- float32 noise -- so only its
  bit-exact digest was regenerated.
- Preserve APB2 root annotation tables and feature relations while annotating or persisting
  ProteoBench results, so independently authored long-form annotations survive the workflow.
- Added `apb-proteobench benchmark` to annotate and score an existing APB2 result.
- Added `apb-proteobench run` for the complete vendor-table → APB2 → FASTA peptide check → ProteoBench workflow. It scores the level named in the module settings as the vendor table reports it.
- Quantitative aggregation is deliberately outside this package. `apb-aggregate` is reached only
  as a separate CLI step, so `apb-proteobench` declares no dependency on it and its only APB
  dependencies are `apb2` and `apb-fasta`.
- Packaged all 11 current ProteoBench module TOMLs with stable catalogue names, checksums, and
  upstream Apache-2.0 provenance. The eight quantitative HYE/HY modules used by legacy APB are
  validated and loadable; plasma, de novo, and entrapment are retained with explicit unsupported
  status for planned integration.

## 0.1.0 — 2026-09-01

- Created the separately released ProteoBench integration for APB2.
- Added strict module-TOML sample annotation with normalized configuration provenance.
- Migrated the legacy golden-tested HYE calculation and added configuration-driven HY support.
- Added explicit diagnostic and scoring protocols, APB result persistence, Cyclopts commands,
  and a verbose Loguru score summary.
- Pinned CI and documentation builds to the compatible APB2 boundary at
  `1d904e664ce1322f852235730287695bfed487b5`.
- Recorded the migrated legacy calculation, golden fixture, and mapper baseline as
  `anndata_bridge/legacy@dd734b1c9d51dfa7e367363d96f8453439fa6235`.
- Aligned local and GitHub Pages documentation builds with the workspace Zensical toolchain.

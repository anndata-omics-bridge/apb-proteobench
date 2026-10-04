# Changes

## 2026-10-04

- `run quant`, `run entrapment` and `benchmark` accept `--scores PATH`, writing the ProteoBench datapoint JSON in the layout of the upstream datapoints, one per q-value kind for entrapment, so APB scores can be compared with ProteoBench's.
- The CLI imports the package only through `api.py`: it moves to `apb_proteobench.cli.app` (console script unchanged) with its presentation, timing and result-performance writers beside it in `cli/`. `ProteoBenchAnalysisResult.submission(layer_name)` returns the submission content the CLI used to assemble from internals; `api.py` exports `LoadedModule` and `SubmissionContent`.
- `api.py` exports `load_module`. `ProteoBenchAnalyzer`, `EntrapmentAnalyzer` and `load_module` drop the `/` and `*` signature markers; every call that worked before still works.
- **Breaking:** apb-proteobench reads FASTA only through apb-fasta and no longer depends on protein-fasta. `EntrapmentAnalyzer(module, fasta)` takes the `FastaAnnotator` of the entrapment FASTA and derives its pairs, so `fasta_pairs` stays internal; the CLI reads FASTA with `FastaAnnotator.read`. The entrapment module drops `mapping_file`, which only located ProteoBench's pair file, from its settings, provenance and packaged TOML, and `run entrapment`'s `load_module` timing now measures module loading. `api.py` exports `load_packaged_module` and `load_packaged_entrapment_module`, which build both analyzers from packaged modules. `available_modules()` and `packaged_module_names()` are gone; the configuration guide's tables list the modules, and the API guide no longer offers `analyze_level()` as a second, lower-level path.

## 2026-10-03

- **Breaking:** `run entrapment` drops `--pairs`: labels and pairs come from the entrapment FASTA, which reproduces ProteoBench's pair file exactly (see [compatibility](docs/compatibility.md)). `fasta_pairs(proteins)` replaces `read_pairs(path)`; `pair` becomes text, the unmodified pair index plus each modification site. The FASTA argument of `run quant` and `run entrapment` may be a protein-fasta database Parquet file.
- Sample annotation runs through `apb2.api.AnnotationCompiler(unmatched="error")` with the module samples as a prolfquapp table keyed by `raw_file`. The annotated obs gains `sample_name` and `condition`, no longer a duplicate `raw_file`; APB2's per-level annotation report moves from the `proteobench` to the `prolfquapp` metadata section, and the module record stays under `proteobench.provenance.annotation`.
- **Breaking:** `ProteoBenchAnalyzer(layers=...)` takes abundance layer names; `None`, the default, scores every abundance layer. `selection=` and `ProteoBenchAnalysisResult.selection` are gone, and scoring provenance records `layers` instead of `selection_mode`. apb2 is imported only through `apb2.api`.
- Use APB2 abundance selections and typed var roles; missing abundance roles fail instead of selecting primary.

## Unreleased

- `entrapment_dia_astral` ships an SDRF for its three HeLa runs; the TOML names it with `sdrf`, and `EntrapmentModuleSettings` records `sdrf` and `sdrf_sha256`, so entrapment provenance identifies it.

- **Breaking:** species and contaminants come from apb-fasta's per-feature FASTA matches (`varm["fasta_validation"]`: `fasta_matching_organisms`, `fasta_matches_contaminant`) instead of the `protein_assignment` column mapped through ProteoBench's `mapper.csv`. Scoring therefore requires peptide verification and fails without it. Deleted `mapper.csv`, `map_reported_proteins` and the level record `protein_mapping`, whose remaining `species_mapper` duplicated the root configuration. The ProteoBot content hash drops the mapper checksum and moves to version 2.

- **Breaking:** `run` is now `apb-proteobench run quant`, with unchanged options; `apb-proteobench run entrapment` scores ProteoBench's entrapment module (`entrapment_dia_astral`) from vendor files, a FASTA and ProteoBench's pair file. `EntrapmentAnalyzer` and `EntrapmentAnalysisResult` are the in-memory API; `load_packaged_entrapment_module` and `read_pairs` load the module and pair file. Each precursor q-value kind that apb-catalog's `proteobench_entrapment` set offers is scored (DIA-NN `q_value`, `library_q_value`, `global_q_value`). Ported from ProteoBench `1147290`; equal q-values share a rank, which differs from ProteoBench only for tied entrapment/target pairs. apb-catalog becomes a dependency.

- Documented the remaining differences from ProteoBench in `docs/compatibility.md`: feature exclusion by the `protein_assignment` column rather than ProteoBench's per-tool `Proteins` column, and the caller-chosen scored layer, with the plasma measurements.

- **Breaking:** `run` and `benchmark` score one abundance layer: `--layer NAME`, default `X` for the APB primary layer. `--x` is removed, and the CLI no longer scores every abundance layer by default; the Python API's `ALL_ABUNDANCE_LAYERS` still does. `--result-performance` therefore needs no extra selector, and exit status 2 for conflicting layer options is gone.

- Every quantitative module's scores now include ProteoBench's plasma metrics, ported from `QuantDatapointPYE` at v0.18.0 with unchanged key names and formulas: spike-in error, per-species and spike-in counts, HUMAN dynamic range and HUMAN error per cutoff, plus four top-level projections. Spike-ins are all species other than HUMAN. Compatibility version moves to 0.18.0 (`cd2a8f00`), whose HYE scoring equals 0.17.0; scoring method `proteobench-compatible` becomes version 2. `ScoreConfig` requires the module `species`. Across 198 corpus intermediates the new metrics equal upstream within 3e-15 relative.

- `dia_plasma` is now a supported module: a TOML plus SDRF for the 12 timsTOF runs of PYE9 (Distler et al. 2025), `max_nr_observed = 12`, the HYE FASTA, and the Custom-format run names as aliases.

- `benchmark` and `run --module` accept a packaged module name such as `dda_qexactive`; a value ending in `.toml` is still read as a module TOML path.

- **Breaking:** a module is now a TOML plus the SDRF-Proteomics table it names (`sdrf = "…"`). The SDRF owns the sample design and species quantities; expected A/B ratios are computed from `characteristics[spiked compound]` quantities. The TOML replaces `species_expected_ratio`/`species_mapper` with `[species.<NAME>]` (`organism`, `suffix`, `color`) and `[[samples]]` with SDRF rows plus `[run_aliases]`. `SampleSettings.raw_file_alias` becomes the list `raw_file_aliases`. Provenance `source` gains `sdrf` (name, checksum) and format `proteobench-module-toml-sdrf`. Resolved settings, matching and scores of all eight packaged modules are unchanged, except `dda_peptidoform`, whose `raw_file` is now the QExactive raw-file name; its former `abundance_*` identifiers remain aliases. sdrf-pipelines joins the dev group to validate the packaged SDRFs.

- **Breaking:** ProteoBot export now names its JSON and dataset with a versioned SHA-256 content hash of the selected abundance layer and scoring inputs. Historical SHA-1 submissions are not matched by the new identity. Scoring no longer computes or persists an upload hash; the optional bundle writer owns it. Species mean and median reductions use Polars while condition reductions remain NumPy.

- **Breaking:** `run` uses one `--software` parameter-grammar hint; removed `--params-software`. APB2 restricts result recognition to the hinted vendor and its declared quantification software (including FragPipe → DIA-NN), or recognizes all vendors when no hint is supplied. Scoring and compatibility exports are unchanged.

- `run --timings-dir DIR` optionally publishes three independent version-1 tool-timing JSON files for in-process APB2 conversion, FASTA verification, and ProteoBench benchmarking. Existing timing targets are refused before conversion; scored results and compatibility exports are unchanged.

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

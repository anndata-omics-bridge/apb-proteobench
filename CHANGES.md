# Changes

## Unreleased

- Added `apb-proteobench convert`, which starts at a vendor table and composes APB2's packaged-rule
  detection and compiler/parser API directly.
- Added `convert_vendor_result()` for single-level h5ad and all-compatible-level h5mu conversion,
  returning the parsed in-memory APB2 result.
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

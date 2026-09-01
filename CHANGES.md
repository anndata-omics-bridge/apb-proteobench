# Changes

## 0.1.0 — 2026-09-01

- Created the separately released ProteoBench integration for APB2.
- Added strict module-TOML sample annotation with normalized configuration provenance.
- Migrated the legacy golden-tested HYE calculation and added configuration-driven HY support.
- Added explicit diagnostic and scoring protocols, APB result persistence, Cyclopts commands,
  and a verbose Loguru score summary.
- Pinned CI and documentation builds to the compatible APB2 boundary at
  `ff496ddee363df01e06c14538d55a60968256a9e`.
- Recorded the migrated legacy calculation, golden fixture, and mapper baseline as
  `anndata_bridge/legacy@dd734b1c9d51dfa7e367363d96f8453439fa6235`.
- Aligned local and GitHub Pages documentation builds with the workspace Zensical toolchain.

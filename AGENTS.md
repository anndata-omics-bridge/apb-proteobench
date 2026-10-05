# APB ProteoBench — agent rules

The closest `AGENTS.md` wins. Explicit user instructions override this file.

## Verified commands

| Task | Command |
| --- | --- |
| Synchronize | `uv sync --group dev` |
| Format | `.venv/bin/ruff format src tests && .venv/bin/ruff check --fix src tests` |
| Lint | `.venv/bin/ruff check src tests` |
| Typecheck | `.venv/bin/pyright` |
| Dependencies | `.venv/bin/deptry .` |
| Tests | `.venv/bin/pytest --cov --cov-branch` |
| Build | `uv build && .venv/bin/twine check dist/*` |
| Full gate | `make check` |

## Architecture

- `api.py` exposes in-memory ProteoBench analysis over canonical `ParsedLevels`; it owns no paths,
  APB2 compilation, FASTA loading, persistence, logging, or CLI behavior.
- `cli/` composes APB2 conversion, FASTA verification, ProteoBench analysis, persistence, and
  presentation through `api.py` only; presentation never reopens a result.
- `integration.py` is the only module translating `ParsedLevels` into calculation inputs or
  attaching calculated outputs.
- `workflow.py` owns the client-side diagnostic and scoring protocols and composes concrete
  calculations explicitly.
- `calculation/` contains no AnnData, MuData, result I/O, logging, or CLI behavior.
- `configuration/` is the inward Pydantic/TOML/SDRF boundary and imports no outer package module; it reads SDRFs only through `apb2.api.SdrfSource`.
- `cli/result_performance.py` and `cli/timings.py` write the CLI's compatibility and timing files
  from completed values; they import neither the CLI application nor its presentation.
- APB2 is an inward dependency through its public facades; APB2 must never import this package.
- Import `ParsedLevels` from `apb2.api` in the public API. Only the CLI may import APB2 compilation,
  quantification-level, and result-I/O operations; only integration may import lower-level result
  types and projection helpers.
- The only permitted APB dependencies are `apb2` and `apb-fasta`. Never import `apb_aggregate` or
  any other sibling APB tool; reach aggregation through the `apb-aggregate` CLI as a separate step.
- Keep HYE and HY configuration-driven. Do not add benchmark-name branches or preset catalogues.
- Preserve the tool-owned layout: root `proteobench.provenance.annotation` and `.scoring` hold common configuration; level `proteobench.annotation` holds matching evidence and `.scoring[quantity_name]` holds each layer's roles, scores and `varm` reference. Scoring provenance is schema 3. No `layers` wrapper, `X` score alias or duplicated selection list; keep reversible layer-name escaping. See [result layout](docs/results.md).

## Compatibility policy

- Pre-1.0 breaking API, CLI, and persisted-schema changes are explicitly allowed when they produce a cleaner public design.
- Do not add compatibility aliases, migration facades, or duplicate persistence layouts unless the user explicitly requests them.
- Record breaking changes in `CHANGES.md`, update callers and documentation in the same change, and fail explicitly on unsupported legacy result layouts.

## Code conventions

- Fully annotate every function and method in `src/` and `tests/`, including
  private functions, callbacks, generators, fixtures, and special methods.
- Standard Pyright strict and Ruff are mandatory. Do not create baselines,
  blanket exclusions, file-wide ignores, or unqualified `# type: ignore`.
- Ruff is the sole formatter and linter. Do not add Black, isort, Flake8, mypy,
  or another overlapping formatter/type checker.
- Keep `__init__.py` empty and import from defining modules inside this package. Other anndata_bridge packages import this one only from `apb_proteobench.api`, and it imports them only from theirs. The CLI imports this package only from `apb_proteobench.api` too.
- Use Google-style docstrings for public APIs and the configured 100-character
  line length.

## Dependency rules

### MUST

- Declare every imported runtime dependency directly in `[project.dependencies]`.
- Put tests, linting, typing, building, and documentation tools in dependency
  groups; optional user-facing capabilities belong in extras.
- Update `pyproject.toml` and run `make check`; the repository commits no `uv.lock`.

### SHOULD

- Prefer the standard library, then an existing direct dependency, then a small,
  maintained, typed dependency.
- Keep source independent of test, build, documentation, and CLI-only packages.

### MUST NOT

- Depend on unpinned branches or undeclared transitive dependencies.
- Add parallel manifests, lockfiles, formatters, type checkers, or test runners.
- Silence a dependency or typing defect instead of fixing its source.

## Workflow

1. Preserve unrelated worktree changes.
2. Add or update focused tests with each behavioral change.
3. Run the smallest relevant check while iterating.
4. Run `make check` before handoff and report its actual result.

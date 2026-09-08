# APB ProteoBench

ProteoBench annotation, mixed-species diagnostics, and scoring for storage-neutral APB2
results. HYE and HY use the same configuration-driven calculation.

Run the complete workflow directly from vendor files:

```bash
apb-proteobench run report.tsv proteins.fasta \
    --params search-parameters.txt \
    --module module_settings.toml \
    --software spectronaut \
    --output results/scored.h5mu
```

This one call converts every compatible level with APB2, verifies modification-stripped peptide sequences against the FASTA database, applies the ProteoBench sample design, calculates diagnostics and scores, and writes one final `results/scored.h5mu`. It scores the level named in `module_settings.toml` as the vendor table reports it and performs no quantitative aggregation.

The equivalent reusable workflow calls the existing tools directly:

```bash
apb2 convert report.tsv --params search-parameters.txt --software spectronaut --output results/all
apb-fasta verify-peptides results/all.h5mu proteins.fasta --output results/fasta-checked.h5mu
apb-aggregate ion protein sum results/fasta-checked.h5mu results/aggregated.h5mu
apb-proteobench benchmark results/aggregated.h5mu module_settings.toml results/scored.h5mu
```

The `apb-aggregate` step is optional and belongs to a separate tool: include it only when the scored level must be derived from a lower one, chaining one call per source level. Omit it and `benchmark` reads `results/fasta-checked.h5mu` directly. APB ProteoBench declares only APB2 and APB FASTA as APB dependencies and reaches aggregation solely as a subprocess. Use the staged form when intermediate results should be inspectable or reusable. `benchmark` combines ProteoBench annotation and scoring for an existing APB2 result. Fine-grained `annotate` and `score` commands remain available. Scoring writes feature diagnostics to `varm["proteobench"]`, and writes scores and provenance to `uns["apb"]["proteobench"]`. See the [documentation](docs/index.md).

The package owns all 11 ProteoBench module TOMLs: the eight quantitative HYE/HY modules used by
legacy APB and the newer plasma, de novo, and entrapment documents. The latter three are explicitly
packaged for planned support but are not accepted by the current quantitative scorer. Stable names,
support status, and the Python loading API are documented under
[module configuration](docs/configuration.md).

## Development

```bash
uv sync --group dev
make check
.venv/bin/pre-commit install --hook-type pre-commit --hook-type pre-push
```

All Python commands run from the synchronized project `.venv`.

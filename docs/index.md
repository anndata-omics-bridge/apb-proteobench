# APB ProteoBench

`apb-proteobench` adds ProteoBench experiment annotation, diagnostics, and scoring to APB2
results without coupling the scientific calculation to AnnData or MuData.

The workflow can start directly from a vendor table:

```bash
apb-proteobench convert report.tsv \
    --params search-parameters.txt \
    --software spectronaut \
    --output converted
apb-proteobench annotate converted.h5mu module_settings.toml annotated.h5mu
apb-proteobench score annotated.h5mu scored.h5mu --verbose
```

`convert` uses APB2's packaged rules and compiler/parser API. With no level it writes every
compatible level to h5mu; with a level such as `ion` it writes one h5ad. The `annotate` and `score`
commands also accept APB2 results produced separately, including h5ad, h5mu, Parquet, and DuckDB.

The annotation command checks that the module describes every observation exactly once. It
stores the validated module configuration with the dataset. The score command therefore needs
no preset name and no second module argument. HYE and HY are two module configurations consumed
by the same calculation.

The resulting feature table lives in `varm["proteobench"]`. Aggregate scores, method identities,
role resolution, compatibility versions, and mapper provenance live in
`uns["apb"]["proteobench"]`.

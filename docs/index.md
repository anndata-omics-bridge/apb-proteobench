# APB ProteoBench

`apb-proteobench` adds ProteoBench experiment annotation, diagnostics, and scoring to APB2
results without coupling the scientific calculation to AnnData or MuData.

The normal workflow is explicit:

```bash
apb2 convert ...
apb-proteobench annotate converted.h5mu module_settings.toml annotated.h5mu
apb-proteobench score annotated.h5mu scored.h5mu --verbose
```

The annotation command checks that the module describes every observation exactly once. It
stores the validated module configuration with the dataset. The score command therefore needs
no preset name and no second module argument. HYE and HY are two module configurations consumed
by the same calculation.

The resulting feature table lives in `varm["proteobench"]`. Aggregate scores, method identities,
role resolution, compatibility versions, and mapper provenance live in
`uns["apb"]["proteobench"]`.

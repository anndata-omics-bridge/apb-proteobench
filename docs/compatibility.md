# Differences from ProteoBench

apb-proteobench reproduces ProteoBench's quantification scores from APB2 results. The calculation is the same; the inputs can differ in two places: which protein column decides feature exclusion, and which abundance layer is scored. Both are listed here with their measured effect.

## Feature exclusion uses the protein-assignment column

Scoring drops contaminant features (a `Cont_` protein) and features of more than one species (`min_count_multispec`). APB decides both from the column the APB2 vendor rule gives the `protein_assignment` role. ProteoBench decides them from each tool's column mapped to `Proteins`.

| Tool | ProteoBench `Proteins` | APB `protein_assignment` |
| --- | --- | --- |
| DIA-NN | `Protein.Ids` | `Protein.Group` |
| FragPipe | `Protein` plus `Mapped Proteins` | `Protein` |
| AlphaDIA | `genes` | `pg.proteins` (2.x) |
| MaxQuant, PEAKS, Sage, Spectronaut, i2MassChroQ, WOMBAT | same column | same column |

DIA-NN's `Protein.Ids` lists every protein containing the peptide; `Protein.Group` is DIA-NN's parsimonious group, which lists several accessions only for proteins the observed precursors cannot distinguish. A precursor whose `Protein.Ids` include a contaminant or a second species is therefore dropped by ProteoBench but kept by APB when parsimony removed that protein from the group.

Measured on one DIA-NN 2.3 plasma submission (22,504 precursors):

- APB keeps 657 precursors that ProteoBench drops; ProteoBench keeps 4 that APB drops
- 458 of the 657 have a `Cont_` protein only in `Protein.Ids`
- 193 have a second species only in `Protein.Ids`
- Over the 12 cutoffs, `nr_feature` is up to 4.3 % higher in APB; spike-in and HUMAN errors differ by up to 1.0 %

Across the 19 DIA-NN and FragPipe plasma submissions, scored on `Precursor_Quantity`, the stored ProteoBench datapoints differ by up to 5.2 % in `nr_feature`, 2.3 % in HUMAN error and 1.3 % in spike-in error; AlphaDIA and PEAKS are identical.

## The scored layer is chosen by the caller

`run` and `benchmark` score one abundance layer: `--layer NAME`, by default `X`, the APB2 rule's primary layer. ProteoBench fixes the quantity column per module and tool, and its plasma module departs from its other DIA modules:

| Tool | ProteoBench plasma | ProteoBench other DIA modules | APB primary (`X`) |
| --- | --- | --- | --- |
| DIA-NN | `Precursor.Quantity` | `Precursor.Normalised` | `Precursor_Normalised` |
| Spectronaut | `FG.Quantity` | `EG.TotalQuantity (Settings)` | `FG_Quantity` |

To reproduce ProteoBench's plasma datapoints from DIA-NN results, pass `--layer Precursor_Quantity`. With `X`, DIA-NN 2.x plasma results differ by more than 50 %, because DIA-NN's normalisation factors differ between the two conditions there.

## Smaller differences

- Two-species (HY) modules write no `nr_quantified_ECOLI`; ProteoBench's plasma code always writes it, as 0
- Per-species keys exist only for the module's species; key names and formulas otherwise follow ProteoBench 0.18.0

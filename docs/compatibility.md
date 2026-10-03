# Differences from ProteoBench

apb-proteobench reproduces ProteoBench's quantification scores from APB2 results. The calculation is the same; the inputs can differ in two places: which protein column decides feature exclusion, and which abundance layer is scored. Both are listed here with their measured effect.

## Feature exclusion uses FASTA matches

Scoring drops contaminant features and features of more than one species (`min_count_multispec`). APB decides both from apb-fasta's `varm["fasta_validation"]`: every FASTA protein containing the feature's peptide gives `fasta_matching_organisms` and `fasta_matches_contaminant`. ProteoBench decides them from each tool's column mapped to `Proteins`.

For DIA-NN, whose `Protein.Ids` lists every protein containing the peptide, both rules agree. Tools whose `Proteins` column is a protein group or a tool-specific list lose features that APB now recognises as shared with a contaminant or a second species.

Measured against the stored ProteoBench datapoints of 219 corpus submissions, `nr_feature` at cutoff 1:

- DIA-NN and FragPipe (DIA-NN quant): median difference 0.1–0.3 %; before this rule, 0.5–3.9 %
- Plasma DIA-NN: 2.9 % before, 0.34 % after
- MaxQuant, AlphaDIA, Spectronaut, PEAKS, Sage, AlphaPept, FragPipe (DDA), quantms: previously within 0.7 % and mostly exact, now up to 3.3 %
- MaxQuant: its 0.2 % of features without a FASTA match are its reversed decoys (`REV__`), which ProteoBench excludes as decoys too
- i2MassChroQ: 21–24 %, because APB2's i2MassChroQ rule keeps modifications in `ProForma_peptide`, so modified peptides match no protein

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

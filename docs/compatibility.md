# Differences from ProteoBench

apb-proteobench reproduces ProteoBench's quantification scores from APB2 results. The calculation is the same; the inputs can differ in two places: which protein column decides feature exclusion, and which abundance layer is scored. Both are listed here with their measured effect.

## Feature exclusion uses FASTA matches

Scoring drops decoys, contaminant features and features of more than one species (`min_count_multispec`). APB decides species and FASTA contaminants from apb-fasta's `varm["fasta_validation"]`: every FASTA protein containing the feature's peptide gives `fasta_matching_organisms` and `fasta_matches_contaminant`. Decoys, and contaminants the software flags or adds itself such as MaxQuant's `CON__`, come from apb2's `apb_Decoy` and `apb_Contaminant`. ProteoBench decides them from each tool's column mapped to `Proteins`.

`--proteins vendor` decides as ProteoBench does, from the proteins the vendor reports: apb2's `fasta_accessions` column, which apb-fasta resolves in the FASTA to `reported_member_organisms` and `reported_members_contaminant`, plus `apb_Contaminant`. ProteoBench reads the same proteins from its `Proteins` column and the species from their `_HUMAN`-style suffix, mapping bare accessions through its `mapper.csv`; the FASTA replaces that mapping. The column differs from ProteoBench's for five rules: AlphaDIA 1.12 and 2.x read `proteins` where ProteoBench reads `genes`, which misses 4,942 and 7,338 `Cont_` precursors; AlphaPept and Spectronaut 19–21 list the same contaminants.

For DIA-NN, whose `Protein.Ids` lists every protein containing the peptide, both rules agree. Tools whose `Proteins` column is a protein group or a tool-specific list lose features that APB now recognises as shared with a contaminant or a second species.

Measured against the stored ProteoBench datapoints of 219 corpus submissions, `nr_feature` at cutoff 1:

- DIA-NN and FragPipe (DIA-NN quant): median difference 0.1–0.3 %; before this rule, 0.5–3.9 %
- Plasma DIA-NN: 2.9 % before, 0.34 % after
- MaxQuant, AlphaDIA, Spectronaut, PEAKS, Sage, AlphaPept, FragPipe (DDA), quantms: previously within 0.7 % and mostly exact, now up to 3.3 %
- MaxQuant: its 0.2 % of features without a FASTA match are its reversed decoys (`REV__`), which ProteoBench excludes as decoys too
- i2MassChroQ: 21–24 %, because APB2's i2MassChroQ rule keeps modifications in `ProForma_peptide`, so modified peptides match no protein

## The scored layer is chosen by the module or the caller

`run` and `benchmark` score one abundance layer: `--layer NAME`, by default `X`, the module's layer. That is the APB2 rule's primary layer, unless the module sets `abundance_normalisation`; then it is the `intensity` layer with that normalisation in apb-catalog's `abundance` set. ProteoBench fixes the quantity column per module and tool, and its plasma module departs from its other DIA modules:

| Tool | ProteoBench plasma | ProteoBench other DIA modules | APB primary (`X`) |
| --- | --- | --- | --- |
| DIA-NN | `Precursor.Quantity` | `Precursor.Normalised` | `Precursor_Normalised` |
| Spectronaut | `FG.Quantity` | `EG.TotalQuantity (Settings)` | `FG_Quantity` |

`dia_plasma` sets `abundance_normalisation = "none"`, so `X` scores DIA-NN's `Precursor_Quantity`, as ProteoBench does. The primary layer would differ by more than 50 % on DIA-NN 2.x plasma results, because DIA-NN's normalisation factors differ between the two conditions there. Where the catalogue names no un-normalised layer, such as for PEAKS, AlphaDIA and Spectronaut, the primary layer is scored and the level summary entry `abundance_normalisation` needs attention.

## Entrapment pairs come from the FASTA

ProteoBench labels each precursor and finds its partner in a separate pair file. APB derives both from the entrapment FASTA, which it already reads to verify peptides:

- Label: protein_fasta's `is_entrapment` flags the `sp|<peptide>_p_target|…` entries
- Unmodified pair: the FASTA lists every target directly before its entrapment, so the pair is the entry position halved
- Modified forms: the pair file pairs each one with the partner form carrying the same modifications on the same occurrence of each residue, e.g. the second Met oxidised; APB's pair key is the unmodified pair plus those sites
- Measured on the 2026 pair file: this key reproduces all 2,573,718 pairs, with no two forms sharing a key
- Guard: scoring stops when two adjacent FASTA entries are not a target followed by a same-length entrapment
- Difference: the pair file lists only the modified forms ProteoBench generated; APB pairs any Oxidation and Carbamidomethyl combination, which the pair-file version refused. None of the 12 entrapment submissions has such a form: all scored under that version.

## Smaller differences

- Two-species (HY) modules write no `nr_quantified_ECOLI`; ProteoBench's plasma code always writes it, as 0
- Per-species keys exist only for the module's species; key names and formulas otherwise follow ProteoBench 0.18.0

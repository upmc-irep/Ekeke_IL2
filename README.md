# Ekeke IL2 dataset and scripts

This repository contains the processed dataset from the paper: "Intrapleural IL-2 Expressing Oncolytic Virotherapy Mediates Superior Antitumor Efficacy and T Cell Receptor Diversity in Malignant Pleural Disease" by Ekeke, et al.

## Sample naming 

A Tumor and Spleen sample was collected from each mouse. Each sample was amplified twice, and the repeat
is indicated by the "_R" in the sample name. In some cases, a reamplification was required: "_Reamp".

The two digit number indicates the "Mouse ID" and S/T indicates Spleen (S) or Tumor (T).

For analysis where repeats were not needed, the amplification with the higher RNA receptor
count (number of UMIs that mapped to a T cell receptor) was used. That list of samples can be found
in [assets/best_samples.csv](assets/best_samples.csv).

## Diversity calculations

Example command used for calculating diversity for all TRB samples.

```bash
./src/diversity_metrics.py $(find csv data/cdr3_pep | grep TRB)
```

## Shared CDR3s between Tumor and Spleen

Command used to calculate sharing percentages.

```bash
./src/spleen_tumor_shared.py assets/grouped_samples.csv
```

## Processed data

The processed data can be found in [processed_data](processed_data). The file format is the "full" preset (default) from mixcr's exportClones command. For more info, see the [mixcr documentation](https://mixcr.readthedocs.io/en/master/export.html). 

To extract all samples, you can run this command:

```bash
for file in data/exported_clones/*.tar.gz; do tar -zxvf $file; done
```

The CDR3 peptide + UMI count files are also included, which were used to calculate CDR3 peptide sharing and diversity.


# Eutherian Corticogenesis

This is a repository for the analysis performed in the paper [Eutherian Enhancers Important in Cortical Development](https://)

## Tissue Extraction

## Published Data

## Cell Ranger Parameters

### Mouse

  --transcriptome=/fast/users/newmana_c/work/10X/refdata-gex-mm10-2020-A \
  --sample=m_mouse \
  --expect-cells=10000 \
  --jobmode=slurm \
  --maxjobs=100 \
  --jobinterval=1000 \
  --include-introns


### Opossum

Opossum transcriptome was generated....

Hi 


  --transcriptome=/fast/users/newmana_c/work/10X/Monodelphis.domestica_genome \
  --expect-cells=10000 \
  --jobmode=slurm \
  --maxjobs=100 \
  --jobinterval=1000 \
  --include-introns

## Analysis in R

### Dependencies

### Normalization & Integration

Integration, normalization and scaling of data output by cellranger runs was integrated using the [opossum_seurat_E14_E16_E18_.Rmd](opossum_seurat_E14_E16_E18_.Rmd) script.

Cells from all samples belong to the pyramidal lineage are then subset and saved in a file 'sc_seurat.PN.sct.rds'.



### Plotting & Differential Expression

### Gene Ontology of Species Specific Differences by Cluster



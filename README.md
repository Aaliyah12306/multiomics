# Multiomics analysis records

This repository records analysis code and verified execution outcomes. Raw datasets and credentials are excluded.

## Initial HPC test

- Job 64293457 failed during the R version check, before reading input.
- Job 64294397 used native R 4.1.1 but stopped because Seurat was unavailable.
- Neither attempt generated figures.
- Next: install compatible packages in a persistent personal R 4.1 library on compute resources, validate, then rerun in a new directory.

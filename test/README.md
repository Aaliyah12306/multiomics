# Native R 4.1.1 test

Jobs 64293457 and 64294397 failed preflight (R version mismatch and missing Seurat respectively). Neither generated figures.

Job 64296101 was submitted with 16 CPUs, 128000 MB memory and a 2-hour limit to install persistent native R 4.1.1 packages, then plot only if package validation succeeds. Installation success is not yet confirmed.

The isolated library is activated with `source ~/test/activate_native_r411.sh` on an allocated compute node. It does not modify the R 4.4.1 library. The installation uses the 2024-06-01 Posit CRAN snapshot and archived Matrix 1.6-5; installed versions are recorded upon success.

`feature_plot.R` is a publication copy: the private input path was replaced with `INPUT_RDS`; this copy itself has not been run. Set that environment variable, then run `Rscript --vanilla feature_plot.R NEW_OUTPUT_DIRECTORY`. The original executed scripts remain on HPC. It plots five genes on the existing harmony.umap; res0.1 metadata is validated, not recomputed. No raw data or credentials are stored here.

Sources: https://docs.hpc.sjtu.edu.cn/app/compilers_and_languages/r.html and R's personal-library mechanism. The HPC manual demonstrates compute-node allocation and Conda; it does not prescribe this native-module package recipe.

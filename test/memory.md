# WholeTissue0810 neuronal markers

Job 64308172 submitted. Results not yet validated.

Input: WholeTissueData_cluster_0810.rds (private HPC location omitted).
Resources: 1 node, 16 CPUs, 128000 MB, 2 hours.
Runtime: official R4.4.1 container, /usr/local/bin/Rscript, existing personal library. Environment and package checks precede data loading.
Analysis: UCHL1, STMN2, PRPH, NDRG4, VIP; existing harmony.umap and unique res0.1 metadata. No reclustering or cell-type filtering.
Expression definition: normalized RNA data > 0 divided by all object cells, with complete coverage required.
Expected outputs: combined five-panel PDF, expression_proportions.csv, validation.txt.
HPC project-relative paths: test/feature_plot.R, test/execution/job.slurm, test/logs, test/results, test/memory.md.
Public script replaces private input/output paths with placeholders; analysis logic unchanged, redacted copy itself not rerun. No raw data or credentials published.
User authorized ongoing updates to these memory records. Next: query existing job and validate results; do not resubmit blindly.

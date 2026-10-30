# Installation follow-up

Job 64296101 failed before plotting because the native R module lacks lattice, a Matrix dependency. No figures were generated.

A new script in a new HPC run directory installs lattice before Matrix, then Seurat and plotting dependencies. Job 64296308 was submitted with the same 16 CPU / 128000 MB / 2-hour limit. The persistent library and activation script remain the designated reusable environment. Package installation and figures are not yet validated. Previous run files are retained.

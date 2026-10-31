# cc-perturbseq-manuscript-analysis

Figure-reproduction curation is maintained on the `curation/figure-reproduction-v1` branch.

## Hippo input recovery and Figures 3–4

After fetching `curation/figure-reproduction-v1` to Sherlock, open
[`notebooks/preprocessing/Rebuild_Hippo_Timecourse.ipynb`](notebooks/preprocessing/Rebuild_Hippo_Timecourse.ipynb)
and use **Run All** in the original analysis environment on a high-memory compute node.
It rebuilds the two missing Hippo H5ADs, validates them, updates local paths, and
runs the affected figure notebooks. [Run instructions and provenance](docs/HIPPO_REBUILD.md).

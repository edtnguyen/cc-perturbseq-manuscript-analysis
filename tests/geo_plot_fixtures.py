from pathlib import Path
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse
from cc_perturbseq import geo_submission as geo

def make_hippo_inputs(root):
    from cc_perturbseq import hippo_projection as hp
    from cc_perturbseq.figure3_panels import FIGURE3F_GENES
    from cc_perturbseq.figure4_panels import FIGURE4C_GENES
    root.mkdir(parents=True)
    genes = list(dict.fromkeys(g for groups in (FIGURE3F_GENES, FIGURE4C_GENES)
                              for group in groups.values() for g in group))
    genes += ["OTHER1", "OTHER2"]
    rng = np.random.default_rng(42)
    conditions = list(hp.COND_PCT_MT_CUTOFFS)
    names = [f"{condition}_rep{rep}" for condition in conditions for rep in (1, 2) for _ in range(3)]
    x = sparse.csr_matrix(rng.poisson(3, size=(len(names), len(genes))).astype(np.float32))
    obs = pd.DataFrame({"sample_name": names, "sublibID": ["s1", "s2", "s8"] * (len(names)//3),
                        "biological_sample": names, "pct_counts_mt": 0.1, "pct_counts_ribo": 20.0,
                        "total_counts": np.asarray(x.sum(axis=1)).ravel(),
                        "n_genes_by_counts": x.getnnz(axis=1)},
                       index=[f"hippo_{i}" for i in range(len(names))])
    rna = ad.AnnData(X=x, obs=obs, var=pd.DataFrame({"gene": genes}, index=[f"ENSG{i}" for i in range(len(genes))]))
    rna.layers["counts"] = x.copy()
    programs = [f"P{i}" for i in range(50)]
    raw = rng.gamma(2, 1, size=(len(names), 50)).astype(np.float32)
    usage = raw / raw.sum(axis=1, keepdims=True)
    p = ad.AnnData(X=usage, obs=hp.add_projection_metadata(obs),
                   var=pd.DataFrame({"program": programs, "program_id_zero_based": np.arange(50)}, index=programs))
    p.obsm[hp.RAW_KEY], p.obsm[hp.USAGE_KEY] = raw, usage
    p.obsm[hp.CLR_KEY] = hp.clr_transform_usage(usage)
    p.uns.update({"cnmf_program_names": programs, "projection_modes": [hp.MODE],
                  "X_projection_mode": hp.MODE, "genes_used_for_projection": genes})
    spectra = root / "reference.tsv"
    pd.DataFrame(rng.gamma(2, 1, size=(50, len(genes))), index=programs, columns=genes).to_csv(spectra, sep="\t")
    p.uns["reference_spectra_sha256"] = geo.file_checksums(spectra)["sha256"]
    qpath, ppath = root / "query.h5ad", root / "projection.h5ad"
    rna.write_h5ad(qpath)
    p.write_h5ad(ppath)
    return qpath, ppath, spectra



def make_hippo_geo(path):
    path = Path(path)
    sources = make_hippo_inputs(path.parent / "hippo_source")
    bundle = geo.build_hippo_bundle(*sources)
    path.parent.mkdir(parents=True, exist_ok=True)
    bundle.write_h5mu(path)
    return sources

"""Real HDF5 roundtrips and original panel producers on small synthetic inputs."""
from pathlib import Path
import json

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite

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


def make_mix_inputs(root):
    from cc_perturbseq.figure_s02_parse_mix import PRIMER_ORDER, summarize_guide_pairs
    expression = root / "splitpipe_out/WT/sub1/all-sample/DGE_filtered"
    expression.mkdir(parents=True)
    rng = np.random.default_rng(73)
    samples = [primer for primer in PRIMER_ORDER for _ in range(12)]
    lineage = [v for _ in PRIMER_ORDER for v in (["K562"]*6 + ["hiPSC"]*6)]
    ids = pd.Index([f"cell_{i:03d}" for i in range(len(samples))], name="cell_id")
    x = sparse.csr_matrix(rng.poisson(2, (len(samples), 60)).astype(np.float32))
    mmwrite(expression / "count_matrix.mtx", x)
    pd.DataFrame({"gene_name": ["MT-CO1"] + [f"GENE{i}" for i in range(59)],
                  "gene_id": [f"ENSG{i}" for i in range(60)]}).to_csv(expression / "all_genes.csv", index=False)
    pd.DataFrame({"bc_wells": ids + "__s1", "sample": samples}).to_csv(expression / "cell_metadata.csv", index=False)
    guide_names = ["K562_g1", "K562_g2", "iPSC_g1", "iPSC_g2", "other_control"]
    reference = root / "splitpipe_out/CRISPR/sub1/all-sample/guide_RNAs_unfiltered/all_guides.csv"
    reference.parent.mkdir(parents=True)
    pd.DataFrame({"gene_name": guide_names, "guide_sequence": ["ACGT"]*5}).to_csv(reference, index=False)
    rows = []
    for index, cell in enumerate(ids):
        if index == len(ids)-1:  # Real zero-guide cell must survive integration.
            continue
        for guide in guide_names:
            rows.append({"cell": cell+"__s1", "gene": "gRNA_"+guide, "guide": guide,
                         "count": int(rng.integers(1, 12))})
    pairs = pd.DataFrame(rows)
    counts = root / "splitpipe_out/CRISPR_combined/grnaxcell_count.tsv"
    counts.parent.mkdir(parents=True)
    pairs[["cell", "gene", "count"]].to_csv(counts, sep="\t", index=False)
    metadata = pd.DataFrame({"sample": samples, "cell_type": lineage}, index=ids)
    frozen = summarize_guide_pairs(pairs[["cell", "guide", "count"]], metadata)
    frozen_path = root / "accepted_S2_cell_metrics.csv"
    frozen.to_csv(frozen_path)
    for primer in PRIMER_ORDER:
        path = root / f"splitpipe_out/WT/sub1/{primer}/report/tscp_counts_subsampled.csv"
        path.parent.mkdir(parents=True)
        path.write_text(",hg38\n10,5\n20,9\n30,12\n40,14\n")
    return frozen_path


def test_hippo_roundtrip_reproduces_original_panels_without_source_files(tmp_path):
    sources = make_hippo_inputs(tmp_path / "private")
    output = tmp_path / geo.HIPPO_FILENAME
    report = geo.export_hippo(*sources, output)
    assert report["status"] == "PASS"
    for path in sources:
        path.unlink()
    bundle = md.read_h5mu(output)
    assert set(bundle.mod) == {"rna", "cnmf"}
    assert bundle.mod["cnmf"].varm["reference_gene_spectra_tpm"].shape[0] == 50
    del bundle
    rendered = tmp_path / "public_figures"
    report = geo.validate_submission(output, render_dir=rendered)
    assert len(report["reproduced_reference_tables"]) == 12
    assert len(list(rendered.glob("*.svg"))) == 4
    assert output.with_suffix(".h5mu.md5").is_file()


def test_mix_roundtrip_sparse_metrics_and_all_panels_without_private_files(tmp_path, monkeypatch):
    import shutil
    from cc_perturbseq import figure_s02_parse_mix_data as data
    root = tmp_path / "private"
    frozen = make_mix_inputs(root)
    output = tmp_path / geo.MIX_FILENAME
    # Freeze accepted labels; an export supplied those labels must not recluster.
    def forbidden(*args, **kwargs):
        raise AssertionError("Private source access or reclustering during public reproduction")
    monkeypatch.setattr(data, "assign_transcriptomic_cell_types", forbidden)
    geo.export_parse_mix(root, output, frozen_metrics_csv=frozen)
    shutil.rmtree(root)
    monkeypatch.setattr(data, "read_expression_and_qc", forbidden)
    report = geo.validate_submission(output, render_dir=tmp_path / "public_figures")
    assert report["status"] == "PASS"
    assert len(report["reproduced_reference_tables"]) == 6
    assert len(list((tmp_path / "public_figures").glob("*.svg"))) == 5
    bundle = md.read_h5mu(output)
    assert sparse.issparse(bundle.mod["crispr"].X)
    assert bundle.mod["crispr"].X[-1].sum() == 0
    assert "cell_047" in bundle.mod["rna"].obs_names


def test_tampered_counts_rejected(tmp_path):
    frozen = make_mix_inputs(tmp_path / "private")
    output = tmp_path / geo.MIX_FILENAME
    geo.export_parse_mix(tmp_path / "private", output, frozen_metrics_csv=frozen)
    bundle = md.read_h5mu(output)
    bundle.mod["crispr"].X.data[0] += 1
    tampered = tmp_path / "tampered.h5mu"
    bundle.write_h5mu(tampered)
    with pytest.raises(ValueError, match="checksum mismatch"):
        geo.validate_submission(tampered)


def test_submission_notebooks_valid_and_compile():
    import nbformat
    root = Path(__file__).resolve().parents[1]
    for name in ("Prepare_GEO_Processed_Data.ipynb", "Test_GEO_Processed_Data.ipynb"):
        nb = nbformat.read(root / "notebooks/submission" / name, as_version=4)
        nbformat.validate(nb)
        for index, cell in enumerate(nb.cells):
            if cell.cell_type == "code":
                compile(cell.source, f"{name}:{index}", "exec")

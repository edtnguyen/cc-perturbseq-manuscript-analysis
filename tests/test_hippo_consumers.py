"""Guard against silent cell loss, skipped QC, and shifted Hippo programs."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scanpy as sc
from scipy import sparse

from cc_perturbseq.figure3_panels import (
    FIGURE3F_D3_CONDITIONS,
    FIGURE3F_GENES,
    _hippo_aligned_metadata,
    _hippo_panel_mask,
    _hippo_program_values,
    prepare_figure3f_dotplots,
)
from cc_perturbseq.figure4_panels import (
    FIGURE4C_D12_CONDITIONS,
    FIGURE4C_GENES,
    prepare_figure4c_dotplots,
)
from cc_perturbseq.hippo_projection import add_projection_metadata


@pytest.fixture
def pair():
    conditions = tuple(FIGURE3F_D3_CONDITIONS) + tuple(FIGURE4C_D12_CONDITIONS)
    genes = sorted({g for groups in (FIGURE3F_GENES, FIGURE4C_GENES)
                    for group in groups.values() for g in group}) + ["OTHER_GENE"]
    rows = []
    for condition in conditions:
        # One retained cell, one excluded sublibrary, and one failed-QC cell.
        for replicate in (1, 2, 3):
            rows.append({"sample_name": f"{condition}_rep{replicate}",
                         "sublibID": "s8" if replicate == 2 else "s1",
                         "pct_counts_mt": 0.1,
                         "pct_counts_ribo": 8.0 if replicate != 3 else 0.0,
                         "biological_sample": f"sample_{len(rows)}"})
    obs = pd.DataFrame(rows, index=[f"cell_{i}" for i in range(len(rows))])
    counts = (np.arange(len(obs) * len(genes)).reshape(len(obs), len(genes)) % 17 + 1).astype(float)
    counts[:, -1] *= 25.0  # Full-library normalisation must include this non-marker.
    obs["total_counts"] = counts.sum(axis=1)
    obs["n_genes_by_counts"] = len(genes)
    query = ad.AnnData(sparse.csr_matrix(counts), obs=obs,
                      var=pd.DataFrame({"gene": genes}, index=[f"ENSG{i}" for i in range(len(genes))]))
    programs = [f"P{i}" for i in range(50)]
    projected = ad.AnnData(np.full((len(obs), 50), 1/50, dtype=np.float32),
                          obs=add_projection_metadata(obs),
                          var=pd.DataFrame({"program": programs,
                                            "program_id_zero_based": np.arange(50)}, index=programs))
    values = np.arange(50, dtype=float) - 24.5
    projected.obsm["X_cnmf_clr_top_specific_200"] = np.tile(values, (len(obs), 1))
    # Reordering is safe only if alignment really uses cell identity.
    projected = projected[::-1].copy()
    return query, projected


def test_full_coverage_required_even_for_cells_later_filtered(pair):
    query, projected = pair
    with pytest.raises(ValueError, match="same cells"):
        _hippo_aligned_metadata(query, projected[1:].copy())


def test_metadata_reconstruction_preserves_condition_and_qc(pair):
    query, projected = pair
    projected.obs = query.obs.loc[projected.obs_names].copy()
    obs = _hippo_aligned_metadata(query, projected)
    assert obs.index.equals(projected.obs_names)
    keep = _hippo_panel_mask(obs, FIGURE3F_D3_CONDITIONS)
    assert keep.sum() == len(FIGURE3F_D3_CONDITIONS)
    assert not obs.loc[keep, "sublibID"].eq("s8").any()
    assert obs.loc[keep, "sample_name"].str.endswith("_rep1").all()


def test_metadata_disagreement_and_false_string_qc_are_rejected(pair):
    query, projected = pair
    projected.obs.iloc[0, projected.obs.columns.get_loc("pct_counts_mt")] = 0.2
    with pytest.raises(ValueError, match="metadata disagree"):
        _hippo_aligned_metadata(query, projected)
    projected.obs["pct_counts_mt"] = query.obs.loc[projected.obs_names, "pct_counts_mt"]
    projected.obs["qc_pass_condition_specific"] = "False"
    with pytest.raises(ValueError, match="QC definition"):
        _hippo_aligned_metadata(query, projected)


def test_required_condition_must_have_retained_cells(pair):
    query, projected = pair
    obs = _hippo_aligned_metadata(query, projected)
    obs.loc[obs["condition"].eq(FIGURE3F_D3_CONDITIONS[0]), "qc_pass_condition_specific"] = False
    with pytest.raises(ValueError, match="no cells"):
        _hippo_panel_mask(obs, FIGURE3F_D3_CONDITIONS)


def test_program_column_metadata_is_required_and_checked(pair):
    _, projected = pair
    observed = _hippo_program_values(projected, ("P2", "P34"))
    np.testing.assert_allclose(observed[0], [-22.5, 9.5])
    projected.var["program_id_zero_based"] = np.arange(1, 51)
    with pytest.raises(ValueError, match="column metadata disagree"):
        _hippo_program_values(projected, ("P2", "P34"))


def test_labelled_program_matrix_uses_labels_instead_of_array_positions(pair):
    _, projected = pair
    key = "X_cnmf_clr_top_specific_200"
    values = pd.DataFrame(projected.obsm[key], index=projected.obs_names,
                          columns=projected.var_names).iloc[:, ::-1]
    projected.obsm[key] = values
    np.testing.assert_allclose(_hippo_program_values(projected, ("P2", "P34"))[0], [-22.5, 9.5])


@pytest.mark.parametrize("builder,conditions", [
    (prepare_figure3f_dotplots, FIGURE3F_D3_CONDITIONS),
    (prepare_figure4c_dotplots, FIGURE4C_D12_CONDITIONS),
])
def test_actual_panel_builder_with_plain_raw_query(pair, tmp_path, monkeypatch, builder, conditions):
    query, projected = pair
    query_path, projected_path = tmp_path / "query.h5ad", tmp_path / "projected.h5ad"
    # Exercise the fallback path absent in the old consumers.
    projected.obs = query.obs.loc[projected.obs_names].copy()
    query.write_h5ad(query_path)
    projected.write_h5ad(projected_path)
    captured = []
    def capture(adata, **kwargs):
        captured.append((adata.copy(), kwargs))
        return object()
    monkeypatch.setattr(sc.pl, "dotplot", capture)
    _, _, gene_source, program_source, metadata = builder(query_path, projected_path)
    assert gene_source["n_cells"].eq(1).all()
    assert program_source["n_cells"].eq(1).all()
    assert metadata["n_common_cells"] == query.n_obs
    assert captured[0][1]["standard_scale"] == "var"
    assert captured[1][1]["layer"] == "clr_top_specific_200"
    chosen = query.obs["sample_name"].str.replace(r"_rep\d+$", "", regex=True).isin(conditions)
    chosen &= query.obs["sample_name"].str.endswith("_rep1")
    totals = query.obs.loc[chosen, "total_counts"]
    target_sum = np.median(totals.to_numpy())
    gene_adata = captured[0][0]
    raw = query[gene_adata.obs_names].X.toarray()
    expected = np.log1p(raw / raw.sum(axis=1, keepdims=True) * target_sum)
    positions = [list(query.var["gene"]).index(gene) for gene in gene_adata.var_names]
    np.testing.assert_allclose(gene_adata.X.toarray(), expected[:, positions], rtol=1e-6)

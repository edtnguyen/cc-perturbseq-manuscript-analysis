"""Bounded numerical parity against the recovered comparison producer."""
from __future__ import annotations

import ast
import json
from pathlib import Path
import tempfile
import unittest
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import non_negative_factorization

from cc_perturbseq import hippo_projection as hp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "recovery/raw_downloads/reproject_cnmf_compare_gene_sets.ipynb"


def original_helpers():
    """Execute only original scientific helpers, not notebook data loads/writes."""
    notebook = json.loads(SOURCE.read_text())
    wanted = {"read_spectra", "build_collapsed_query_matrix", "project_fixed_spectra", "clr_transform_usage"}
    tree = ast.parse("".join(notebook["cells"][5]["source"]))
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in wanted]
    assert len(nodes) == len(wanted)
    namespace = {"np": np, "pd": pd, "sparse": sparse,
                 "non_negative_factorization": non_negative_factorization,
                 "TARGET_SUM": 1e6, "SKLEARN_MAX_ITER": 1000,
                 "SKLEARN_TOL": 1e-5, "RANDOM_STATE": 0}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"), namespace)
    return namespace, notebook


class ProjectionParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original, cls.source = original_helpers()

    def test_specificity_ranking_union_order_and_missing_gene_selection(self):
        rng = np.random.default_rng(914)
        # Repeated low integer weights deliberately create specificity/weight ties.
        ref = pd.DataFrame(rng.integers(0, 4, size=(3, 240)),
                           index=["P0", "P1", "P2"],
                           columns=[f"G{i:03}" for i in rng.permutation(240)])
        query = pd.Index(list(ref.columns[7:]) + [ref.columns[10]])
        tree = ast.parse("".join(self.source["cells"][9]["source"]))
        nodes = []
        collecting = False
        for node in tree.body:
            names = {target.id for target in node.targets if isinstance(target, ast.Name)} if isinstance(node, ast.Assign) else set()
            if "gene_total" in names:
                collecting = True
            if collecting:
                nodes.append(node)
            if "top_specific_union_genes" in names:
                break
        namespace = {"np": np, "pd": pd, "W_ref": ref,
                     "query_genes_set": set(query), "TOP_SPECIFIC_N": 200}
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"), namespace)
        self.assertEqual(hp.top_specific_genes(ref, query), namespace["top_specific_union_genes"])

    def test_duplicate_symbols_full_library_totals_and_float32_match_source(self):
        counts = sparse.csr_matrix(np.array([[1.1, 2.2, 10.0, 0, 4], [4, 1, 3.4, 2, 8]], dtype=np.float64))
        ids = pd.Index(["B", "A", "A", "C", "outside"])
        selected = ["C", "A", "absent", "B"]
        totals = np.array([101.25, 202.5], dtype=np.float64)
        expected, expected_genes = self.original["build_collapsed_query_matrix"](counts, ids, selected, totals)
        actual, genes = hp.build_collapsed_query_matrix(counts, ids, selected, totals)
        self.assertEqual(genes, expected_genes)
        np.testing.assert_array_equal(actual.toarray(), expected.toarray())
        self.assertEqual(actual.dtype, np.float32)
        self.assertAlmostEqual(float(actual[0, 1]), np.float32(12.2 / 101.25 * 1e6), places=2)

    def test_fixed_reference_numerical_solution_matches_source(self):
        rng = np.random.default_rng(816)
        reference = rng.uniform(0.1, 3.0, size=(4, 22))
        counts = sparse.csr_matrix(rng.uniform(0.01, 20, size=(12, 22)).astype(np.float32))
        original_ref = reference.copy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            expected_raw, expected_norm, expected_iter = self.original["project_fixed_spectra"](counts, reference.copy(), "test")
            raw, norm, diagnostics = hp.project_fixed_spectra(counts, reference)
        np.testing.assert_array_equal(reference, original_ref)
        np.testing.assert_array_equal(raw, expected_raw)
        np.testing.assert_array_equal(norm, expected_norm)
        self.assertEqual(diagnostics["n_iter"], expected_iter)
        np.testing.assert_array_equal(hp.clr_transform_usage(norm), self.original["clr_transform_usage"](expected_norm))

    def test_clr_floor_handles_zero_entries_without_delta_subtraction(self):
        usages = np.array([[0, 0.3, 0.7], [0.2, 0.3, 0.5]], dtype=np.float32)
        np.testing.assert_array_equal(hp.clr_transform_usage(usages), self.original["clr_transform_usage"](usages))
        self.assertTrue(np.isfinite(hp.clr_transform_usage(usages)).all())
        with self.assertRaises(ValueError):
            hp.clr_transform_usage(np.zeros((2, 3)))

    def test_original_qc_boundaries_are_annotation_only(self):
        obs = pd.DataFrame({"sample_name": ["D1_MPS_rep1"] * 4,
                            "pct_counts_mt": [2.99, 3.0, 2.0, 2.0],
                            "pct_counts_ribo": [2.31, 3.0, 2.3, 2.29]},
                           index=["c0", "c1", "c2", "c3"])
        result = hp.add_projection_metadata(obs)
        self.assertEqual(len(result), 4)
        self.assertEqual(result["qc_pass_condition_specific"].tolist(), [True, False, False, False])
        self.assertEqual(result["condition"].tolist(), ["D1_MPS"] * 4)
        self.assertNotIn("condition", obs)

    def test_one_based_reference_rows_map_to_P0(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectra.tsv"
            frame = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=[1, 2], columns=["A", "B", "C"])
            frame.to_csv(path, sep="\t")
            pd.testing.assert_frame_equal(hp.read_spectra(path, strict_historical=False), self.original["read_spectra"](path))


class ProjectionArtifactTests(unittest.TestCase):
    def test_roundtrip_provenance_resume_and_full_cell_correspondence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            query_path, reference_path, output = root / "query.h5ad", root / "spectra.tsv", root / "projected.h5ad"
            rng = np.random.default_rng(2026)
            obs = pd.DataFrame({
                "sample_name": ["D1_MPS_rep1"] * 6 + ["D2_LM+TDI_rep2"] * 6,
                "pct_counts_mt": [1.0] * 12, "pct_counts_ribo": [4.0] * 12,
                "total_counts": np.arange(12) + 1000.5, "sublibID": ["s1"] * 12,
                "biological_sample": ["original"] * 12,
            }, index=[f"cell{i}_s1" for i in range(12)])
            var = pd.DataFrame({"gene": [f"G{i}" for i in range(15)] + ["G1"]},
                               index=[f"ENS{i}" for i in range(16)])
            query = ad.AnnData(sparse.csr_matrix(rng.uniform(0, 20, size=(12, 16))), obs=obs, var=var)
            query.write_h5ad(query_path)
            pd.DataFrame(rng.uniform(0.01, 3, size=(4, 15)), index=[1, 2, 3, 4],
                         columns=[f"G{i}" for i in range(15)]).to_csv(reference_path, sep="\t")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = hp.project_top_specific_200(query_path, reference_path, output, strict_historical=False)
            self.assertFalse(result["reused"])
            self.assertTrue(result["validation"]["exact_cell_order"])
            self.assertEqual(result["validation"]["projected_shape"], [12, 4])
            self.assertEqual(result["validation"]["qc_passing_cells"], 12)
            resumed = hp.project_top_specific_200(query_path, reference_path, output, strict_historical=False)
            self.assertTrue(resumed["reused"])
            projected = ad.read_h5ad(output)
            subset_path = root / "missing_cell.h5ad"
            projected[:-1].copy().write_h5ad(subset_path)
            with self.assertRaisesRegex(ValueError, "every query cell"):
                hp.validate_projected(query_path, subset_path, strict_historical=False)
            reversed_path = root / "reversed.h5ad"
            projected[::-1].copy().write_h5ad(reversed_path)
            with self.assertRaisesRegex(ValueError, "every query cell"):
                hp.validate_projected(query_path, reversed_path, strict_historical=False)
            projected.obsm[hp.CLR_KEY][0, 0] += 1
            tampered = root / "wrong_clr.h5ad"
            projected.write_h5ad(tampered)
            with self.assertRaisesRegex(ValueError, "original CLR"):
                hp.validate_projected(query_path, tampered, strict_historical=False)
            with self.assertRaises(FileExistsError):
                hp.project_top_specific_200(query_path, reference_path, tampered, strict_historical=False)


if __name__ == "__main__":
    unittest.main()

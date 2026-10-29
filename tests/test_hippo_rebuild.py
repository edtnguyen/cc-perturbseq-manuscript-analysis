"""Bounded regression tests for the production Hippo expression producer."""
from __future__ import annotations

import os
from pathlib import Path
import tempfile
import unittest

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from cc_perturbseq.hippo_rebuild import (
    HISTORICAL_SAMPLE_COUNTS,
    INPUT_RELATIVE_TEMPLATE,
    SAMPLE_MAPPING,
    combined_receipt_path,
    preflight_sublibraries,
    rebuild_combined,
    validate_combined,
)


def write_source_fixture(root: Path) -> list[ad.AnnData]:
    """Eight small sources; each boundary row tests an explicit production rule."""
    sources = []
    for i in range(1, 9):
        # Fractional raw quantification counts must remain unchanged.
        totals = np.array([2000.25, 1300.75, 2000, 2000, 1200, 50000, 2000, 2000, 0.0])
        values = np.column_stack([totals * 0.25, totals * 0.5, totals * 0.25])
        obs = pd.DataFrame({
            "n_genes_by_counts": [1301, 1302, 1200, 12000, 1301, 1301, 1301, 1301, 0],
            "total_counts": totals,
            "pct_counts_mt": [5, 4, 5, 5, 5, 5, 10, 5, np.nan],
            "pct_counts_ribo": [2.5] * 8 + [np.nan],
            "biological_sample": pd.Categorical(["1-1", " 1-2 "] + ["1-1"] * 5 + ["empty", "empty"]),
            "plate": ["plate1"] * 9,
            "well": ["A1"] * 9,
            "pool": ["Hippo"] * 9,
        }, index=pd.Index([f"cell{j}" for j in range(9)], name="barcode"))
        var = pd.DataFrame({"gene": ["CDH5", "ETV2", "TAGLN"], "mt": [False] * 3},
                           index=pd.Index(["g1", "g2", "g3"], name="gene_id"))
        current = ad.AnnData(sparse.csr_matrix(values), obs=obs, var=var)
        for name, factor in (("ambiguous", 0.1), ("mature", 0.6), ("nascent", 0.3)):
            current.layers[name] = sparse.csr_matrix(values * factor)
        current.obsm["guide_assignment"] = sparse.csr_matrix(np.ones((9, 2)))
        current.obsm["guide_counts"] = sparse.csr_matrix(np.arange(18).reshape(9, 2))
        current.uns["guide_names"] = np.array(["guide1", "guide2"])
        # An outer concat must align by feature identity, including a different order.
        if i == 8:
            current = current[:, [2, 0, 1]].copy()
        path = root / INPUT_RELATIVE_TEMPLATE.format(i=i)
        path.parent.mkdir(parents=True)
        current.write_h5ad(path)
        sources.append(current)
    return sources


class HippoRebuildTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.inputs = self.root / "inputs"
        self.output = self.root / "output" / "Hippo_timecourse_combined.h5ad"
        self.sources = write_source_fixture(self.inputs)

    def tearDown(self):
        self.directory.cleanup()

    def rebuild(self, **kwargs):
        return rebuild_combined(self.inputs, self.output, strict_historical=False, **kwargs)

    def test_boundaries_identity_layers_and_recorded_dimensions(self):
        original_signatures = preflight_sublibraries(self.inputs)
        receipt = self.rebuild()
        result = ad.read_h5ad(self.output)
        self.assertEqual(result.shape, (16, 3))
        self.assertEqual(receipt["pre_qc_shape"], [72, 3])
        self.assertEqual(receipt["post_qc_shape"], [16, 3])
        self.assertFalse(receipt["condition_specific_qc_applied"])
        expected_ids = [f"cell{cell}_s{i}" for i in range(1, 9) for cell in (0, 1)]
        self.assertEqual(result.obs_names.tolist(), expected_ids)
        self.assertEqual(result.obs["sample_name"].astype(str).tolist(), ["D0_hPSC_rep1", "D0_hPSC_rep2"] * 8)
        self.assertEqual(result.obs["biological_sample"].astype(str).tolist(), ["1-1", "1-2"] * 8)
        for i, source in enumerate(self.sources, start=1):
            actual = result[[f"cell0_s{i}", f"cell1_s{i}"], source.var_names]
            np.testing.assert_array_equal(actual.X.toarray(), source.X[:2].toarray())
            for name in ("ambiguous", "mature", "nascent"):
                np.testing.assert_array_equal(actual.layers[name].toarray(), source.layers[name][:2].toarray())
            for name in ("guide_counts", "guide_assignment"):
                np.testing.assert_array_equal(actual.obsm[name].toarray(), source.obsm[name][:2].toarray())
        np.testing.assert_array_equal(result.uns["guide_names"], ["guide1", "guide2"])
        self.assertEqual(original_signatures, preflight_sublibraries(self.inputs))
        self.assertEqual(receipt["validation"]["sample_counts"], {"1-1": 8, "1-2": 8})
        self.assertEqual(self.rebuild()["status"], "reused")

    def test_missing_input_fails_before_any_output(self):
        missing = self.inputs / INPUT_RELATIVE_TEMPLATE.format(i=4)
        missing.unlink()
        with self.assertRaisesRegex(FileNotFoundError, "sub4"):
            self.rebuild()
        self.assertFalse(self.output.parent.exists())

    def test_dataframe_guide_annotations_are_preserved_and_validated(self):
        for i in range(1, 9):
            path = self.inputs / INPUT_RELATIVE_TEMPLATE.format(i=i)
            source = ad.read_h5ad(path)
            for name in ("guide_assignment", "guide_counts"):
                source.obsm[name] = pd.DataFrame(
                    source.obsm[name].toarray(), index=source.obs_names,
                    columns=["guide1", "guide2"],
                )
            source.write_h5ad(path)
        receipt = self.rebuild()
        self.assertEqual(receipt["validation"]["shape"], [16, 3])
        result = ad.read_h5ad(self.output)
        for name in ("guide_assignment", "guide_counts"):
            self.assertIsInstance(result.obsm[name], pd.DataFrame)
            self.assertEqual(result.obsm[name].shape, (16, 2))
            self.assertEqual(result.obsm[name].columns.tolist(), ["guide1", "guide2"])
            self.assertEqual(result.obsm[name].index.tolist(), result.obs_names.tolist())

    def test_unknown_existing_output_is_not_overwritten(self):
        self.output.parent.mkdir()
        self.output.write_bytes(b"existing user file")
        with self.assertRaises(FileExistsError):
            self.rebuild()
        self.assertEqual(self.output.read_bytes(), b"existing user file")

    def test_changed_source_signature_is_not_silently_reused(self):
        self.rebuild()
        changed = self.inputs / INPUT_RELATIVE_TEMPLATE.format(i=2)
        stat = changed.stat()
        os.utime(changed, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
        with self.assertRaises(FileExistsError):
            self.rebuild()

    def test_strict_production_rejects_wrong_cohort_without_writing(self):
        with self.assertRaisesRegex(ValueError, "Pre-QC shape"):
            rebuild_combined(self.inputs, self.output)
        self.assertFalse(self.output.exists())

    def test_normalized_or_corrupt_expression_is_rejected(self):
        self.rebuild()
        with h5py.File(self.output, "r+") as handle:
            handle["X"]["data"][:] /= 2
        with self.assertRaisesRegex(ValueError, "X row sums"):
            validate_combined(self.output, strict_historical=False)
        with h5py.File(self.output, "r+") as handle:
            handle["X"]["data"][0] = np.nan
        with self.assertRaisesRegex(ValueError, "non-finite raw counts"):
            validate_combined(self.output, strict_historical=False)

    def test_unknown_sample_preserves_source_order_of_operations(self):
        path = self.inputs / INPUT_RELATIVE_TEMPLATE.format(i=1)
        source = ad.read_h5ad(path)
        source.obs["biological_sample"] = source.obs["biological_sample"].astype(str)
        source.obs.loc["cell0", "biological_sample"] = "empty "
        source.write_h5ad(path)
        # Source cell 14 does not strip before excluding 'empty'; cell 23 then
        # reveals that the retained label has no sample mapping.
        with self.assertRaisesRegex(ValueError, "Unmapped biological_sample"):
            self.rebuild()
        self.assertFalse(self.output.exists())

    def test_invalid_input_does_not_replace_a_completed_result(self):
        self.rebuild()
        previous = self.output.read_bytes()
        path = self.inputs / INPUT_RELATIVE_TEMPLATE.format(i=1)
        source = ad.read_h5ad(path)
        source.layers["mature"].data[0] = -1
        source.write_h5ad(path)
        with self.assertRaisesRegex(ValueError, "negative or non-finite"):
            self.rebuild(overwrite=True)
        self.assertEqual(previous, self.output.read_bytes())
        self.assertFalse(list(self.output.parent.glob("*.tmp.h5ad")))

    def test_recorded_historical_counts_cover_the_exact_sample_map(self):
        self.assertEqual(sum(HISTORICAL_SAMPLE_COUNTS.values()), 212891)
        self.assertEqual(set(HISTORICAL_SAMPLE_COUNTS), set(SAMPLE_MAPPING))
        self.assertEqual(len(SAMPLE_MAPPING), 54)


if __name__ == "__main__":
    unittest.main()

# Public reproduction environment and GEO upload boundaries

See `GEO_PROCESSED_EXPORTS.md` for the two H5MU schemas and export procedure.

## Existing Sherlock environment

Use the environment that produced the accepted panels. Check that `mudata` is installed before running the exporter. Do not replace the production environment with the public test environment merely to export the data: the existing S2 classification is version-sensitive. Supplying the complete accepted-run per-cell metrics CSV retains its lineage labels without reclustering. Both H5MU objects record the actual export package versions.

## New environment for public-file testing

A Python 3.11 environment with pinned core scientific dependencies is provided in `config/geo-reproduction-requirements.txt`. These versions match the successful synthetic HDF5 roundtrip environment. This is a core-version specification, not a complete transitive lockfile.

From a fresh Python 3.11 virtual environment in the repository:

```bash
python -m pip install -r config/geo-reproduction-requirements.txt
python -m ipykernel install --user --name ccperturbseq-geo --display-name "CC-Perturb-seq GEO"
PYTHONPATH=src python -m pytest -q tests/test_geo_submission.py
```

Open `notebooks/submission/Test_GEO_Processed_Data.ipynb` using that kernel. Set `CCPERTURBSEQ_GEO_DIR` to the local directory containing the two downloaded H5MU objects, or place them in `<repository>/geo_release/`. It uses no Sherlock configuration or external data files. Its temporary Hippo H5ADs are extracted from the public H5MU locally and then deleted.

The automatic tests use synthetic data and must not be confused with validation of the study files. For study validation, run the testing notebook on the actual exported objects and inspect its panel previews.

## Upload only the actual submission files

The two H5MU files are the processed data. Submit the original FASTQs and the required study/sample metadata separately according to GEO's instructions. Do not gzip the H5MU files again: they already use internal HDF5 compression.

The `.md5`, `.sha256`, validation JSON, temporary extracts, rendered test figures, and code ZIP are local audit products, not files to upload automatically. GEO currently requests that optional MD5 values be entered in the metadata spreadsheet's MD5 Checksums tab rather than submitted as separate checksum files. Upload only files listed in the submission metadata; do not upload a ZIP archive of the workflow.

Single-cell GEO sample records normally represent a pool/library with its unique raw-file set, not individual cells. Make the sample/library identifiers in the submission metadata traceable to the cell-level annotations and explicitly describe the relationship between expression libraries, guide libraries, sublibraries and pooled biological conditions. The notebooks do not infer missing experimental metadata or fabricate sample/FASTQ mappings.

Official guidance (checked 30 October 2026):
https://www.ncbi.nlm.nih.gov/geo/info/seq.html

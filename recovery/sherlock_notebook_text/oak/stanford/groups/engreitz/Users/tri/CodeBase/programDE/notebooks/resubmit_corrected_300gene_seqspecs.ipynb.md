# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/resubmit_corrected_300gene_seqspecs.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/resubmit_corrected_300gene_seqspecs.ipynb`
- SHA-256: `c0eab01e6df488ad6f87cd4969d4a45ff2b45aa98ca738ae81774dd41312b031`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Regenerate and replace the 48 corrected 300-gene CC-Perturb-seq seqspec files

This notebook performs a **full structural-template correction** for the 24 GEX
and 24 gRNA seqspec ConfigurationFiles that were originally submitted in May
2026.

## Authoritative correction policy

For each modality, all structural and onlist fields come from the supplied
corrected example:

- GEX: `IGVFFI6906STSO.yaml.gz`
- gRNA: `IGVFFI4006RNQZ.yaml.gz`

For every original sample-specific seqspec, the notebook preserves:

- original `date`;
- original `assay_id`;
- both original `read_id` values;
- every original FASTQ `!File` entry, including accession, filename, size, URL,
  and MD5;
- existing ConfigurationFile accession and alias;
- existing `file_set` and `seqspec_of` links.

The notebook regenerates all 48 YAML and deterministic YAML.GZ files under:

```text
/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/
seqspec_two_template_exact/corrected_templates/
```

It then updates the **existing** 48 ConfigurationFile records in place. It does
not POST duplicate records.

## Safety model

- All source files must match their original submitted MD5 values.
- All corrected files must match the authoritative template after replacing
  only the explicitly preserved sample-specific fields.
- Portal alias, file set, and `seqspec_of` must remain unchanged.
- Unexpected Portal MD5 values or released/archived/revoked records block the
  workflow.
- Metadata dry runs are required before production.
- Production checks use hashes of the generated files and plans.
- Each S3 upload is verified immediately with `head-object` when possible.

## Markdown cell 1

## 1. Configuration

## Code cell 2

```python
from __future__ import annotations

import base64
import copy
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import time
import zipfile
from getpass import getpass
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pandas as pd
import requests
import yaml

# Default Sherlock locations. Environment variables can override these only for
# testing; the defaults are the intended production paths.
SOURCE_ROOT = Path(
    os.environ.get(
        "SEQSPEC_SOURCE_ROOT",
        "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/"
        "seqspec_two_template_exact",
    )
)

OUTPUT_ROOT = Path(
    os.environ.get(
        "SEQSPEC_CORRECTED_ROOT",
        str(SOURCE_ROOT / "corrected_templates"),
    )
)

YAML_DIR = OUTPUT_ROOT / "yaml"
YAML_GZ_DIR = OUTPUT_ROOT / "yaml_gz"
REFERENCE_DIR = OUTPUT_ROOT / "reference_templates"
BACKUP_DIR = OUTPUT_ROOT / "backup_original_files"
PORTAL_BACKUP_DIR = OUTPUT_ROOT / "backup_portal_records"
REPORT_DIR = OUTPUT_ROOT / "reports"
LOG_DIR = OUTPUT_ROOT / "logs"
WORK_DIR = OUTPUT_ROOT / "_work"

for directory in (
    OUTPUT_ROOT,
    YAML_DIR,
    YAML_GZ_DIR,
    REFERENCE_DIR,
    BACKUP_DIR,
    PORTAL_BACKUP_DIR,
    REPORT_DIR,
    LOG_DIR,
    WORK_DIR,
):
    directory.mkdir(parents=True, exist_ok=True)

API_BASE = "https://api.data.igvf.org"
REQUEST_TIMEOUT_SECONDS = 120

IGVF_PYTHON = Path(
    "/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python"
)

SEQSPEC_EXECUTABLE_CANDIDATES = [
    Path(
        "/oak/stanford/groups/engreitz/Users/tri/envs/seqspec-0.4/bin/seqspec"
    ),
    Path(
        "/oak/stanford/groups/engreitz/Users/tri/envs/seqspec/bin/seqspec"
    ),
]

RUN_METADATA_DRY_RUNS = True

PRODUCTION_CONFIRMATION = ""
REQUIRED_PRODUCTION_CONFIRMATION = (
    "REPLACE_48_300GENE_SEQSPEC_FILES_IN_PLACE"
)

# When True, production attempts a Portal download after upload. S3 head-object
# verification remains the immediate authoritative upload check.
TRY_PORTAL_DOWNLOAD_VERIFICATION = False

# Embedded production submission map parsed from the original successful log.
SUBMISSION_RECORDS = json.loads("[{\"i\":\"1\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI9366VIOI\",\"file_set\":\"\\/measurement-sets\\/IGVFDS6519LIUB\\/\",\"md5sum\":\"1732aaccba8653587a41f0752e55f512\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI0829JETG,IGVFFI0873SMHV\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D0_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"2\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI5296LCGA\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS5155REOE\\/\",\"md5sum\":\"8ad8a723246603c4f5c7bb6053c02419\",\"file_size\":\"1549\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI8551DLZI,IGVFFI9572IFKJ\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D0_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"3\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI3664MJLW\",\"file_set\":\"\\/measurement-sets\\/IGVFDS7946PPWP\\/\",\"md5sum\":\"86c7f4571e6301af80ea7cde10c60d75\",\"file_size\":\"1418\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI1612DYPE,IGVFFI8493FUGF\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D0_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"4\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI6014NOEH\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS9937DZRX\\/\",\"md5sum\":\"954679771bcaa7e2c39fb48b5ffbf746\",\"file_size\":\"1547\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI0467PSNK,IGVFFI2082AZES\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D0_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"5\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI9935VMFV\",\"file_set\":\"\\/measurement-sets\\/IGVFDS0186YGHF\\/\",\"md5sum\":\"69a2e1c8aadc5fc8d4173146de6f7e9b\",\"file_size\":\"1414\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI6695QNYC,IGVFFI9053NTKZ\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D1_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"6\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI6690AXDH\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS4763NXDE\\/\",\"md5sum\":\"a462480e775c4eb1e8dbc51f63addd45\",\"file_size\":\"1548\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI3402IDBG,IGVFFI3106VMTS\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D1_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"7\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI2772BUCG\",\"file_set\":\"\\/measurement-sets\\/IGVFDS1346LIZG\\/\",\"md5sum\":\"1e4c8c11b5e0a0e43ab5cb55c0bdaf3c\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI5558MDRQ,IGVFFI8244CRTY\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D1_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"8\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI0810QEPN\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS4792JAEL\\/\",\"md5sum\":\"78f57be06976c5901d196ad202ef962d\",\"file_size\":\"1548\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI3532ZCVD,IGVFFI9691ZIPG\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D1_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"9\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI8048DGMK\",\"file_set\":\"\\/measurement-sets\\/IGVFDS3062PAWU\\/\",\"md5sum\":\"426974f55aba03d3e2f58f6baad6ac3f\",\"file_size\":\"1416\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI3535RMGA,IGVFFI9133VCWO\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D2_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"10\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI8873GOQM\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS9204CZCV\\/\",\"md5sum\":\"0b4d08704f9720c63df43726233711f3\",\"file_size\":\"1549\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI2247UTAX,IGVFFI5633MUBM\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D2_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"11\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI4312QTKI\",\"file_set\":\"\\/measurement-sets\\/IGVFDS3687IPQU\\/\",\"md5sum\":\"0f252bf0c96811a57653199292cdbd79\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI9459ZZMQ,IGVFFI4270HSCI\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D2_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"12\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI9788HFUZ\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS2592NCPO\\/\",\"md5sum\":\"295c189806a4f4d99a58c26ed068fa04\",\"file_size\":\"1548\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI9145BQRD,IGVFFI2507CKOT\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D2_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"13\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI7854RUAJ\",\"file_set\":\"\\/measurement-sets\\/IGVFDS6661ZZZO\\/\",\"md5sum\":\"75e149f184203d27b6e69f8c3b974f51\",\"file_size\":\"1416\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI4758OTQB,IGVFFI6464CXQX\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D3_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"14\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI9066UFCN\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS0523AHQL\\/\",\"md5sum\":\"586c86f0c716f6971fad30ae5c6b3cda\",\"file_size\":\"1548\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI2633HLYP,IGVFFI8658UHCI\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D3_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"15\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI2187SJIK\",\"file_set\":\"\\/measurement-sets\\/IGVFDS7648OGLO\\/\",\"md5sum\":\"504cfb99a4e53e1b6e10fb51f5f63439\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI1136PBFN,IGVFFI3187HIDR\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D3_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"16\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI8158UWTN\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS0313DTTA\\/\",\"md5sum\":\"f31c629f25e92e6fd16978e13f02d580\",\"file_size\":\"1546\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI3766KXGG,IGVFFI9603VRZJ\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D3_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"17\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI4622KMQE\",\"file_set\":\"\\/measurement-sets\\/IGVFDS0576WBRT\\/\",\"md5sum\":\"8010aa911e3700af327315f0c29cd3f6\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI1685UGHT,IGVFFI4812NIJP\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D0_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"18\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI8076FSQO\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS6717HJOD\\/\",\"md5sum\":\"0743cc132b0c42b016f656e9f1ed6302\",\"file_size\":\"1551\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI1509ELYW,IGVFFI4095EXNM\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D0_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"19\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI2175EAGA\",\"file_set\":\"\\/measurement-sets\\/IGVFDS3972WQLN\\/\",\"md5sum\":\"3f1887d0ef53b32117a5b83730bcf003\",\"file_size\":\"1418\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI2803KCCY,IGVFFI8343RQJO\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D0_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"20\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI2434RZUD\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS5213QJBY\\/\",\"md5sum\":\"57a0705a5daee08e617bca2709adbec7\",\"file_size\":\"1551\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI9495FYTZ,IGVFFI0841HSYI\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D0_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"21\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI8770RPJT\",\"file_set\":\"\\/measurement-sets\\/IGVFDS1003MCWV\\/\",\"md5sum\":\"2c5f7f42f522363b9455e7c2aa830653\",\"file_size\":\"1416\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI3189ORTZ,IGVFFI0206VVDK\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D1_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"22\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI1488HGCB\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS4718CSYF\\/\",\"md5sum\":\"73a056d1f142b6222d67adb322c8f365\",\"file_size\":\"1549\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI7823OLFM,IGVFFI6204QQIO\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D1_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"23\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI3154UJVP\",\"file_set\":\"\\/measurement-sets\\/IGVFDS3023PHPP\\/\",\"md5sum\":\"6fb418da684953915bec5872ce77d613\",\"file_size\":\"1414\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI0050MVZE,IGVFFI2411QTXC\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D1_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"24\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI2812VJVV\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS1487DXOE\\/\",\"md5sum\":\"33605d7ca6fccca62e68d858a68f0042\",\"file_size\":\"1550\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI9198TJDQ,IGVFFI6658QMWQ\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D1_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"25\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI2347AMLW\",\"file_set\":\"\\/measurement-sets\\/IGVFDS0918LDGC\\/\",\"md5sum\":\"84707f5bc3c53958ee62b4515e3ec8d8\",\"file_size\":\"1415\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI0665JNTA,IGVFFI0999QMPT\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D2_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"26\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI7718VJVV\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS1745LTCI\\/\",\"md5sum\":\"8065e829ba1f7607f24eb24f6ed0cdb6\",\"file_size\":\"1548\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI4016NGIZ,IGVFFI2998SBIB\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D2_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"27\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI5545EVXX\",\"file_set\":\"\\/measurement-sets\\/IGVFDS4438MNHW\\/\",\"md5sum\":\"a620eab20a00d45ac2b5f8acfd2b196c\",\"file_size\":\"1414\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI2484ZZPZ,IGVFFI9124BSXO\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D2_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"28\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI0654IPLV\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS4706FMSW\\/\",\"md5sum\":\"8b501d304b699945a32ad0bf058bd62b\",\"file_size\":\"1549\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI8861UHPT,IGVFFI7405SFKA\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D2_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"29\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI8139DPYA\",\"file_set\":\"\\/measurement-sets\\/IGVFDS5755WPHZ\\/\",\"md5sum\":\"be53cb2203a9ec5e99b6eae5d50c6fcf\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI2965GQIJ,IGVFFI1664QYJP\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D3_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"30\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI1456NTTN\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS6690DHIT\\/\",\"md5sum\":\"eb9ac220be5183c562452995e03232f2\",\"file_size\":\"1548\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI7830VVBQ,IGVFFI8521CNTM\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D3_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"31\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI8337KFXH\",\"file_set\":\"\\/measurement-sets\\/IGVFDS4408YWND\\/\",\"md5sum\":\"0873189f042bbebc67109cc036358298\",\"file_size\":\"1418\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI9685TYQG,IGVFFI1590FZYI\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D3_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"32\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI6865SJSV\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS8999NJNY\\/\",\"md5sum\":\"bafa1e59ff12ef948edbc436588b22bc\",\"file_size\":\"1549\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI9284YBSZ,IGVFFI8667VZSA\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D3_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"33\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI5741OOFJ\",\"file_set\":\"\\/measurement-sets\\/IGVFDS0124WAPG\\/\",\"md5sum\":\"040016c01790799c41536e418d114c76\",\"file_size\":\"1415\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI0770LMLA,IGVFFI9913JDMS\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D0_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"34\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI3023CUFB\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS0090ZFOF\\/\",\"md5sum\":\"5e14b4005a133c32a89d5524ce15776e\",\"file_size\":\"1547\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI5792YSTV,IGVFFI7768XINN\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D0_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"35\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI9717GCXG\",\"file_set\":\"\\/measurement-sets\\/IGVFDS9673XGGF\\/\",\"md5sum\":\"1fc0e5735f1cdd12dc6980d0487a112f\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI0292PKWS,IGVFFI3456DMRR\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D0_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"36\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI5837FVFN\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS6554CSYJ\\/\",\"md5sum\":\"b20f9df58e5ee5c8bb4eb174cf6364e2\",\"file_size\":\"1548\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI3227OROC,IGVFFI7228ZGTG\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D0_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"37\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI6943RGAA\",\"file_set\":\"\\/measurement-sets\\/IGVFDS3909AASU\\/\",\"md5sum\":\"0f62f6e45f86460b1c317ab4f6bce6ab\",\"file_size\":\"1419\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI8066JZTB,IGVFFI9818VWXH\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D1_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"38\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI6783SCTB\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS0983NLQN\\/\",\"md5sum\":\"3fce6114f140473b1b64f78ac1815ef8\",\"file_size\":\"1549\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI2003RXVG,IGVFFI6693ASNH\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D1_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"39\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI3958VUUO\",\"file_set\":\"\\/measurement-sets\\/IGVFDS3709FLMW\\/\",\"md5sum\":\"8341195801076e2b2c2284af9b2aaeba\",\"file_size\":\"1418\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI4958SRSE,IGVFFI9062DNGV\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D1_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"40\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI4268ZJMH\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS6291BGWR\\/\",\"md5sum\":\"fe36c7811ae1c76c7f7097a682bbb873\",\"file_size\":\"1550\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI7449DXPM,IGVFFI9750FRHU\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D1_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"41\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI4966LKOU\",\"file_set\":\"\\/measurement-sets\\/IGVFDS3333KHQB\\/\",\"md5sum\":\"e7e419696a4d80ccb357b78c9727ee75\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI6246BWHT,IGVFFI1862FOTZ\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D2_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"42\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI9383WWAH\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS1114RWFQ\\/\",\"md5sum\":\"df9da071b5e2b2e9b53c30f97606db9f\",\"file_size\":\"1550\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI5495IAYW,IGVFFI0411AIRZ\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D2_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"43\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI5241GSQF\",\"file_set\":\"\\/measurement-sets\\/IGVFDS8080YPLU\\/\",\"md5sum\":\"ddcd63b44f793d1e24b9fbec6d6e04fe\",\"file_size\":\"1419\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI9466RZXY,IGVFFI2864XLLF\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D2_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"44\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI8680AUTL\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS4732AAFZ\\/\",\"md5sum\":\"b3a095eb3314f27150b4abe8e29e30db\",\"file_size\":\"1548\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI1203DQMX,IGVFFI9532XJBL\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D2_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"45\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep1_GEX_seqspec\",\"accession\":\"IGVFFI3799GGVF\",\"file_set\":\"\\/measurement-sets\\/IGVFDS2398BXCP\\/\",\"md5sum\":\"ea5db335e95a3ab73c14f3cc7a201d1e\",\"file_size\":\"1417\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep1_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI7587JKFC,IGVFFI9686HBHI\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D3_biorep1; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"46\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep1_gRNA_seqspec\",\"accession\":\"IGVFFI8059GUMK\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS0959OQHG\\/\",\"md5sum\":\"58cac223def7545f750a9c3b6220736d\",\"file_size\":\"1549\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep1_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI8543DRYQ,IGVFFI7114OHFF\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D3_biorep1; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"},{\"i\":\"47\",\"modality\":\"GEX\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep2_GEX_seqspec\",\"accession\":\"IGVFFI5429IGZN\",\"file_set\":\"\\/measurement-sets\\/IGVFDS4676VINP\\/\",\"md5sum\":\"14dec7650b72bca549a9458daf3272f3\",\"file_size\":\"1416\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep2_GEX_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI7691VNOQ,IGVFFI6400NHMG\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D3_biorep2; generated using IGVFFI8161HYYS.yaml.gz as template.\"},{\"i\":\"48\",\"modality\":\"gRNA\",\"alias\":\"jesse-engreitz:configuration_file_WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep2_gRNA_seqspec\",\"accession\":\"IGVFFI2335DNRR\",\"file_set\":\"\\/auxiliary-sets\\/IGVFDS8359WXUO\\/\",\"md5sum\":\"a005075ebb1c428ed607256547fc18b4\",\"file_size\":\"1550\",\"submitted_file_name\":\"\\/oak\\/stanford\\/groups\\/engreitz\\/Users\\/tri\\/300genes_IGVF_demux\\/seqspec_two_template_exact\\/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep2_gRNA_seqspec.yaml.gz\",\"seqspec_of\":\"IGVFFI5552ZPHT,IGVFFI7298UQGU\",\"description\":\"WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D3_biorep2; generated using CC-perturb-seq_crispr_IGVFFI9348HZUK.yaml as template.\"}]")

# Embedded corrected reference examples. These are written to Sherlock by the
# notebook and are not required as separate input files.
CORRECTED_REFERENCE_TEXT = {
    "GEX": "!Assay\nseqspec_version: 0.3.0\nassay_id: CC-Peturb-seq\nname: CC-Perturb-seq using 48 barcode set for Round 1 with single illumina multiplex index\ndoi: https://www.protocols.io/view/evercode-wt-mega-v2-2-1-d4rb8v2n\ndate: 08 September 2026\ndescription: split-pool ligation-based for single cell transcriptome and Circular Capture for gRNA \nmodalities:\n- rna\nlib_struct: https://igvf.github.io/seqspec/\nlibrary_protocol: Custom\nlibrary_kit: Custom kit \nsequence_protocol: Illumina NovaSeq X (NTR:0000765)\nsequence_kit: NovaSeq X Series 10B Reagent Kit\nsequence_spec:\n- !Read\n  read_id: IGVFDS3636ULPH_read_1\n  name: Read 1 \n  modality: rna\n  primer_id: truseq_read1\n  min_len: 100\n  max_len: 100\n  strand: pos\n  files:\n  - !File\n    file_id: IGVFFI0289BNWF\n    filename: IGVFFI0289BNWF.fastq.gz\n    filetype: fastq\n    filesize: 28288342565\n    url: https://api.data.igvf.org/sequence-files/IGVFFI0289BNWF/@@download/IGVFFI0289BNWF.fastq.gz\n    urltype: https\n    md5: a2a477e26ace354f89d7623f097dec33\n  - !File\n    file_id: IGVFFI5436FYZN\n    filename: IGVFFI5436FYZN.fastq.gz\n    filetype: fastq\n    filesize: 28628994535\n    url: https://api.data.igvf.org/sequence-files/IGVFFI5436FYZN/@@download/IGVFFI5436FYZN.fastq.gz\n    urltype: https\n    md5: 5ed4c872a0b617ea9a26988594517b9e\n- !Read\n  read_id: IGVFDS3636ULPH_read_2\n  name: Read 2 \n  modality: rna\n  primer_id: truseq_read2\n  min_len: 86\n  max_len: 86\n  strand: neg\n  files:\n  - !File\n    file_id: IGVFFI7081GGDU\n    filename: IGVFFI7081GGDU.fastq.gz\n    filetype: fastq\n    filesize: 37035025159\n    url: https://api.data.igvf.org/sequence-files/IGVFFI7081GGDU/@@download/IGVFFI7081GGDU.fastq.gz\n    urltype: https\n    md5: 7131e16edcb77f673c9d005d2b58677a\n  - !File\n    file_id: IGVFFI9134JTSD\n    filename: IGVFFI9134JTSD.fastq.gz\n    filetype: fastq\n    filesize: 37776346872\n    url: https://api.data.igvf.org/sequence-files/IGVFFI9134JTSD/@@download/IGVFFI9134JTSD.fastq.gz\n    urltype: https\n    md5: 8791f430a281d87bd5d67be7760be6eb\nlibrary_spec:\n- !Region\n  region_id: rna\n  region_type: named\n  name: rna\n  sequence_type: joined\n  sequence: AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATCTXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAAAAAAAAAAAAANNNNNNNNCCACAGTCTCAAGCACGTGGATNNNNNNNNAGTCGTACGCCGATGCGAAACATCGGCCACNNNNNNNNXXXXXXXXXXAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG\n  min_len: 320\n  max_len: 320\n  onlist: null\n  parent_id: null\n  regions:\n  - !Region\n    region_id: illumina_p5\n    region_type: illumina_p5\n    name: Illumina P5\n    sequence_type: fixed\n    sequence: AATGATACGGCGACCACCGAGATCTACAC\n    min_len: 29\n    max_len: 29\n    onlist: null\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: truseq_read1\n    region_type: truseq_read1\n    name: Truseq Single Index Read 1\n    sequence_type: fixed\n    sequence: ACACTCTTTCCCTACACGACGCTCTTCCGATCT\n    min_len: 33\n    max_len: 33\n    onlist: null\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: cdna\n    region_type: cdna\n    name: Read 1 cDNA sequence\n    sequence_type: random\n    sequence: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n    min_len: 100\n    max_len: 100\n    onlist: null\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: poly_A\n    region_type: poly_A\n    name: poly_A\n    sequence_type: fixed\n    sequence: AAAAAAAAAAAAAA\n    min_len: 14\n    max_len: 14\n    onlist: null\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: barcode-1\n    region_type: barcode\n    name: barcode-1\n    sequence_type: onlist\n    sequence: NNNNNNNN\n    min_len: 8\n    max_len: 8\n    onlist: !Onlist\n      file_id: IGVFFI9984SVPZ\n      filename: IGVFFI9984SVPZ.tsv.gz\n      filetype: tsv\n      filesize: 397\n      location: remote\n      url: https://api.data.igvf.org/tabular-files/IGVFFI9984SVPZ/@@download/IGVFFI9984SVPZ.tsv.gz\n      urltype: https\n      md5: 08b4623964ad7a3fea0d70b39e5459eb\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: linker-1\n    region_type: linker\n    name: linker-1\n    sequence_type: fixed\n    sequence: CCACAGTCTCAAGCACGTGGAT\n    min_len: 22\n    max_len: 22\n    onlist: null\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: barcode-2\n    region_type: barcode\n    name: barcode-2\n    sequence_type: onlist\n    sequence: NNNNNNNN\n    min_len: 8\n    max_len: 8\n    onlist: !Onlist\n      file_id: IGVFFI1138MCVX\n      filename: IGVFFI1138MCVX.tsv.gz\n      filetype: tsv\n      filesize: 341\n      location: remote\n      url: https://api.data.igvf.org/tabular-files/IGVFFI1138MCVX/@@download/IGVFFI1138MCVX.tsv.gz\n      urltype: https\n      md5: 7d310eca0f1030f13fec3395d044142d\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: linker-2\n    region_type: linker\n    name: linker-2\n    sequence_type: fixed\n    sequence: AGTCGTACGCCGATGCGAAACATCGGCCAC\n    min_len: 30\n    max_len: 30\n    onlist: null\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: barcode-3\n    region_type: barcode\n    name: barcode-3\n    sequence_type: onlist\n    sequence: NNNNNNNN\n    min_len: 8\n    max_len: 8\n    onlist: !Onlist\n      file_id: IGVFFI1138MCVX\n      filename: IGVFFI1138MCVX.tsv.gz\n      filetype: tsv\n      filesize: 341\n      location: remote\n      url: https://api.data.igvf.org/tabular-files/IGVFFI1138MCVX/@@download/IGVFFI1138MCVX.tsv.gz\n      urltype: https\n      md5: 7d310eca0f1030f13fec3395d044142d\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: umi\n    region_type: umi\n    name: umi\n    sequence_type: random\n    sequence: XXXXXXXXXX\n    min_len: 10\n    max_len: 10\n    onlist: null\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: truseq_read2\n    region_type: truseq_read2\n    name: Truseq Read 2\n    sequence_type: fixed\n    sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n    min_len: 34\n    max_len: 34\n    onlist: null\n    regions: null\n    parent_id: rna\n  - !Region\n    region_id: illumina_p7\n    region_type: illumina_p7\n    name: Illumina P7\n    sequence_type: fixed\n    sequence: ATCTCGTATGCCGTCTTCTGCTTG\n    min_len: 24\n    max_len: 24\n    onlist: null\n    regions: null\n    parent_id: rna\n",
    "gRNA": "!Assay\nseqspec_version: 0.3.0\nassay_id: CC-Peturb-seq\nname: CC-Perturb-seq using 48 barcode set for Round 1 with single illumina multiplex index\ndoi: https://www.protocols.io/view/evercode-wt-mega-v2-2-1-d4rb8v2n\ndate: 08 September 2026\ndescription: split-pool ligation-based for single cell transcriptome and circular capture for guides\nmodalities:\n- crispr\nlib_struct: https://igvf.github.io/seqspec/\nlibrary_protocol: Custom\nlibrary_kit: CC-Perturb-seq\nsequence_protocol: Illumina NovaSeq X (NTR:0000765)\nsequence_kit: NovaSeq X Series 10B Reagent Kit\nsequence_spec:\n- !Read\n  read_id: IGVFDS3967WXNG_read_1\n  name: Read 1\n  modality: crispr\n  primer_id: truseq_read1\n  min_len: 90\n  max_len: 90\n  strand: pos\n  files:\n  - !File\n    file_id: IGVFFI2834DTMA\n    filename: IGVFFI2834DTMA.fastq.gz\n    filetype: fastq\n    filesize: 3088804957\n    url: https://api.data.igvf.org/sequence-files/IGVFFI2834DTMA/@@download/IGVFFI2834DTMA.fastq.gz\n    urltype: https\n    md5: 229a8a21ad365f6001776647b95c7e43\n  - !File\n    file_id: IGVFFI3741PVML\n    filename: IGVFFI3741PVML.fastq.gz\n    filetype: fastq\n    filesize: 3031813408\n    url: https://api.data.igvf.org/sequence-files/IGVFFI3741PVML/@@download/IGVFFI3741PVML.fastq.gz\n    urltype: https\n    md5: 344cf0682d34f6173c242111cf7c7a1f\n- !Read\n  read_id: IGVFDS3967WXNG_read_2\n  name: Read 2\n  modality: crispr\n  primer_id: truseq_read2\n  min_len: 86\n  max_len: 86\n  strand: neg\n  files:\n  - !File\n    file_id: IGVFFI7207CKCB\n    filename: IGVFFI7207CKCB.fastq.gz\n    filetype: fastq\n    filesize: 4542598854\n    url: https://api.data.igvf.org/sequence-files/IGVFFI7207CKCB/@@download/IGVFFI7207CKCB.fastq.gz\n    urltype: https\n    md5: fa9f95fe3f5215086caa5873ffa7ce4e\n  - !File\n    file_id: IGVFFI8642HYJQ\n    filename: IGVFFI8642HYJQ.fastq.gz\n    filetype: fastq\n    filesize: 4458233230\n    url: https://api.data.igvf.org/sequence-files/IGVFFI8642HYJQ/@@download/IGVFFI8642HYJQ.fastq.gz\n    urltype: https\n    md5: 3bf4aed5f8fe97432b9b99cce0ba2837\nlibrary_spec:\n- !Region\n  parent_id: null\n  region_id: crispr\n  region_type: crispr\n  name: crispr\n  sequence_type: joined\n  sequence: AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATCTXXXXXXXAAGCAGTGGTATCAACGCAGAGTGAATGGGNNNNNNNNNNNNNNNNNNNNNGTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGNNNNNNNNCCACAGTCTCAAGCACGTGGATNNNNNNNNAGTCGTACGCCGATGCGAAACATCGGCCACNNNNNNNNXXXXXXXXXXAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG\n  min_len: 296\n  max_len: 297\n  onlist: null\n  regions:\n  - !Region\n    parent_id: crispr\n    region_id: illumina_p5\n    region_type: illumina_p5\n    name: Illumina P5\n    sequence_type: fixed\n    sequence: AATGATACGGCGACCACCGAGATCTACAC\n    min_len: 29\n    max_len: 29\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: truseq_read1\n    region_type: truseq_read1\n    name: Truseq Single Index Read 1\n    sequence_type: fixed\n    sequence: ACACTCTTTCCCTACACGACGCTCTTCCGATCT\n    min_len: 33\n    max_len: 33\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: random_sequence\n    region_type: linker\n    name: random_sequence\n    sequence_type: random\n    sequence: XXXXXXX\n    min_len: 7\n    max_len: 7\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: prefix\n    region_type: custom_primer\n    name: prefix\n    sequence_type: fixed\n    sequence: AAGCAGTGGTATCAACGCAGAGTGAATGGG\n    min_len: 30\n    max_len: 30\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: sgrna_target\n    region_type: gdna\n    name: sgrna_target\n    sequence_type: onlist\n    sequence: NNNNNNNNNNNNNNNNNNNNN\n    min_len: 20\n    max_len: 21\n    onlist: !Onlist\n      file_id: IGVFFI5765HMZH\n      filename: IGVFFI5765HMZH.tsv.gz\n      filetype: tsv\n      filesize: 11046\n      location: remote\n      url: https://api.data.igvf.org/tabular-files/IGVFFI5765HMZH/@@download/IGVFFI5765HMZH.tsv.gz\n      urltype: https\n      md5: cdfa2c279554780920a63641114db587\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: suffix\n    region_type: linker\n    name: suffix\n    sequence_type: fixed\n    sequence: GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAG\n    min_len: 33\n    max_len: 33\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: barcode-1\n    region_type: barcode\n    name: barcode-1\n    sequence_type: onlist\n    sequence: NNNNNNNN\n    min_len: 8\n    max_len: 8\n    onlist: !Onlist\n      file_id: IGVFFI9984SVPZ\n      filename: IGVFFI9984SVPZ.tsv.gz\n      filetype: tsv\n      filesize: 397\n      location: remote\n      url: https://api.data.igvf.org/tabular-files/IGVFFI9984SVPZ/@@download/IGVFFI9984SVPZ.tsv.gz\n      urltype: https\n      md5: 08b4623964ad7a3fea0d70b39e5459eb\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: linker-1\n    region_type: linker\n    name: linker-1\n    sequence_type: fixed\n    sequence: CCACAGTCTCAAGCACGTGGAT\n    min_len: 22\n    max_len: 22\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: barcode-2\n    region_type: barcode\n    name: barcode-2\n    sequence_type: onlist\n    sequence: NNNNNNNN\n    min_len: 8\n    max_len: 8\n    onlist: !Onlist\n      file_id: IGVFFI1138MCVX\n      filename: IGVFFI1138MCVX.tsv.gz\n      filetype: tsv\n      filesize: 341\n      location: remote\n      url: https://api.data.igvf.org/tabular-files/IGVFFI1138MCVX/@@download/IGVFFI1138MCVX.tsv.gz\n      urltype: https\n      md5: 7d310eca0f1030f13fec3395d044142d\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: linker-2\n    region_type: linker\n    name: linker-2\n    sequence_type: fixed\n    sequence: AGTCGTACGCCGATGCGAAACATCGGCCAC\n    min_len: 30\n    max_len: 30\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: barcode-3\n    region_type: barcode\n    name: barcode-3\n    sequence_type: onlist\n    sequence: NNNNNNNN\n    min_len: 8\n    max_len: 8\n    onlist: !Onlist\n      file_id: IGVFFI1138MCVX\n      filename: IGVFFI1138MCVX.tsv.gz\n      filetype: tsv\n      filesize: 341\n      location: remote\n      url: https://api.data.igvf.org/tabular-files/IGVFFI1138MCVX/@@download/IGVFFI1138MCVX.tsv.gz\n      urltype: https\n      md5: 7d310eca0f1030f13fec3395d044142d\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: umi\n    region_type: umi\n    name: umi\n    sequence_type: random\n    sequence: XXXXXXXXXX\n    min_len: 10\n    max_len: 10\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: truseq_read2\n    region_type: truseq_read2\n    name: Truseq Read 2\n    sequence_type: fixed\n    sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n    min_len: 34\n    max_len: 34\n    onlist: null\n    regions: null\n  - !Region\n    parent_id: crispr\n    region_id: illumina_p7\n    region_type: illumina_p7\n    name: Illumina P7\n    sequence_type: fixed\n    sequence: ATCTCGTATGCCGTCTTCTGCTTG\n    min_len: 24\n    max_len: 24\n    onlist: null\n    regions: null\n",
}

print("Source root:", SOURCE_ROOT)
print("Corrected output root:", OUTPUT_ROOT)
print("Embedded submission rows:", len(SUBMISSION_RECORDS))
```

### Text output (stream)

```text
Source root: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact
Corrected output root: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates
Embedded submission rows: 48

```

## Markdown cell 3

## 2. Tagged-YAML, hashing, gzip, diff, and Portal helpers

## Code cell 4

```python
class TaggedMap(dict):
    """Dictionary retaining an explicit seqspec YAML tag."""

    def __init__(self, *args, tag: str | None = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.tag = tag


class TaggedLoader(yaml.SafeLoader):
    pass


def _construct_tagged(loader, tag_suffix, node):
    mapping = loader.construct_mapping(node, deep=True)
    return TaggedMap(mapping, tag="!" + tag_suffix)


TaggedLoader.add_multi_constructor("!", _construct_tagged)


class TaggedDumper(yaml.SafeDumper):
    pass


def _represent_tagged(dumper, data):
    return dumper.represent_mapping(
        data.tag,
        data.items(),
    )


TaggedDumper.add_representer(TaggedMap, _represent_tagged)


def load_tagged_yaml_text(text: str) -> TaggedMap:
    loaded = yaml.load(text, Loader=TaggedLoader)
    if not isinstance(loaded, TaggedMap) or loaded.tag != "!Assay":
        raise RuntimeError("Seqspec root must be an !Assay mapping")
    return loaded


def load_tagged_yaml_gz(path: Path) -> TaggedMap:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return load_tagged_yaml_text(handle.read())


def dump_tagged_yaml(value: TaggedMap) -> str:
    text = yaml.dump(
        value,
        Dumper=TaggedDumper,
        sort_keys=False,
        allow_unicode=True,
        width=120,
        default_flow_style=False,
    )
    if not text.endswith("\n"):
        text += "\n"
    return text


def write_deterministic_gzip(
    source_bytes: bytes,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as raw:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw,
            compresslevel=9,
            mtime=0,
        ) as compressed:
            compressed.write(source_bytes)


def md5_bytes(value: bytes) -> str:
    return hashlib.md5(value).hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def md5_file(path: Path) -> str:
    return md5_bytes(path.read_bytes())


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
    )


def safe_stem(value: object) -> str:
    return re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        str(value),
    ).strip("_")


def read_name(read: dict) -> str:
    return str(read.get("name", "")).strip().lower()


def all_file_ids(spec: dict) -> list[str]:
    output: list[str] = []
    for read in spec.get("sequence_spec", []):
        for file_object in read.get("files", []) or []:
            file_id = str(file_object.get("file_id", "")).strip()
            if file_id:
                output.append(file_id)
    return output


def file_set_accession_from_read_ids(spec: dict) -> str:
    prefixes = set()
    for read in spec.get("sequence_spec", []):
        read_id = str(read.get("read_id", ""))
        match = re.fullmatch(r"(IGVFDS[A-Z0-9]+)_read_[12]", read_id)
        if not match:
            raise RuntimeError(f"Unexpected read_id: {read_id!r}")
        prefixes.add(match.group(1))
    if len(prefixes) != 1:
        raise RuntimeError(f"Read IDs refer to multiple file sets: {prefixes}")
    return next(iter(prefixes))


def normalized_template_object(spec: TaggedMap) -> TaggedMap:
    """Remove only fields explicitly preserved from each old sample."""
    normalized = copy.deepcopy(spec)
    normalized["date"] = "<PRESERVED_ORIGINAL_DATE>"
    normalized["assay_id"] = "<PRESERVED_ORIGINAL_ASSAY_ID>"
    for read in normalized.get("sequence_spec", []):
        read["read_id"] = f"<PRESERVED_{read_name(read)}_ID>"
        read["files"] = ["<PRESERVED_SAMPLE_FASTQ_FILES>"]
    return normalized


def build_corrected_spec(
    old_spec: TaggedMap,
    corrected_template: TaggedMap,
) -> TaggedMap:
    corrected = copy.deepcopy(corrected_template)

    corrected["date"] = old_spec["date"]
    corrected["assay_id"] = old_spec["assay_id"]

    old_reads = {
        read_name(read): read
        for read in old_spec.get("sequence_spec", [])
    }
    corrected_reads = {
        read_name(read): read
        for read in corrected.get("sequence_spec", [])
    }

    if set(old_reads) != {"read 1", "read 2"}:
        raise RuntimeError(f"Old seqspec has unexpected reads: {sorted(old_reads)}")
    if set(corrected_reads) != {"read 1", "read 2"}:
        raise RuntimeError(
            f"Corrected template has unexpected reads: {sorted(corrected_reads)}"
        )

    for key in ("read 1", "read 2"):
        corrected_reads[key]["read_id"] = old_reads[key]["read_id"]
        corrected_reads[key]["files"] = copy.deepcopy(old_reads[key]["files"])

    return corrected


def recursive_diff(old: Any, new: Any, path: str = "$") -> list[dict]:
    rows: list[dict] = []

    if isinstance(old, dict) and isinstance(new, dict):
        keys = list(dict.fromkeys([*old.keys(), *new.keys()]))
        for key in keys:
            child = f"{path}.{key}"
            if key not in old:
                rows.append({"path": child, "old": "<MISSING>", "new": repr(new[key])})
            elif key not in new:
                rows.append({"path": child, "old": repr(old[key]), "new": "<MISSING>"})
            else:
                rows.extend(recursive_diff(old[key], new[key], child))
        return rows

    if isinstance(old, list) and isinstance(new, list):
        maximum = max(len(old), len(new))
        for index in range(maximum):
            child = f"{path}[{index}]"
            if index >= len(old):
                rows.append({"path": child, "old": "<MISSING>", "new": repr(new[index])})
            elif index >= len(new):
                rows.append({"path": child, "old": repr(old[index]), "new": "<MISSING>"})
            else:
                rows.extend(recursive_diff(old[index], new[index], child))
        return rows

    if old != new:
        rows.append({"path": path, "old": repr(old), "new": repr(new)})

    return rows


def region_by_id(spec: dict) -> dict[str, dict]:
    roots = spec.get("library_spec", [])
    if len(roots) != 1:
        raise RuntimeError("Expected exactly one root library region")
    regions = roots[0].get("regions", []) or []
    return {str(region["region_id"]): region for region in regions}


def assert_joined_region_consistency(spec: dict) -> None:
    roots = spec.get("library_spec", [])
    if len(roots) != 1:
        raise RuntimeError("Expected exactly one root library region")
    root = roots[0]
    regions = root.get("regions", []) or []
    joined = "".join(str(region.get("sequence", "")) for region in regions)
    min_len = sum(int(region.get("min_len", 0)) for region in regions)
    max_len = sum(int(region.get("max_len", 0)) for region in regions)

    if joined != root.get("sequence"):
        raise RuntimeError("Root joined sequence does not equal child concatenation")
    if min_len != int(root.get("min_len", -1)):
        raise RuntimeError("Root min_len does not equal child min_len sum")
    if max_len != int(root.get("max_len", -1)):
        raise RuntimeError("Root max_len does not equal child max_len sum")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtnXx", "TGCANtgcanXx"))[::-1]


def expected_read2_pattern(spec: dict) -> str:
    regions = region_by_id(spec)
    ordered = [
        regions["umi"]["sequence"],
        regions["barcode-3"]["sequence"],
        reverse_complement(regions["linker-2"]["sequence"]),
        regions["barcode-2"]["sequence"],
        reverse_complement(regions["linker-1"]["sequence"]),
        regions["barcode-1"]["sequence"],
    ]
    return "".join(map(str, ordered))


def portal_record_url(identifier: str) -> str:
    token = str(identifier).strip()
    if token.startswith("/"):
        return API_BASE.rstrip("/") + token
    return API_BASE.rstrip("/") + "/" + token.strip("/") + "/"


def enter_igvf_credentials() -> None:
    api_key = <REDACTED>"IGVF_API_KEY", "").strip()
    secret_key = os.environ.get("IGVF_SECRET_KEY", "").strip()

    if not api_key:
        <REDACTED> = getpass("IGVF API key: ").strip()
    if not secret_key:
        secret_key = getpass("IGVF secret key: ").strip()

    if not api_key or not secret_key:
        raise RuntimeError("Both IGVF credentials are required")

    os.environ["IGVF_API_KEY"] = api_key
    os.environ["IGVF_SECRET_KEY"] = secret_key
    print("IGVF credentials are set for the current kernel.")


def make_portal_session() -> requests.Session:
    api_key = <REDACTED>"IGVF_API_KEY", "").strip()
    secret_key = os.environ.get("IGVF_SECRET_KEY", "").strip()
    if not api_key or not secret_key:
        raise RuntimeError("Run enter_igvf_credentials() first")

    session = requests.Session()
    session.auth = (api_key, secret_key)
    session.headers.update({"Accept": "application/json"})
    return session


def get_portal_record(
    session: requests.Session,
    identifier: str,
    *,
    frame: str = "object",
    database: bool = True,
) -> dict:
    params = {"format": "json", "frame": frame}
    if database:
        params["datastore"] = "database"
    response = session.get(
        portal_record_url(identifier),
        params=params,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )
    if not response.ok:
        raise RuntimeError(
            f"GET failed for {identifier}: HTTP {response.status_code}\n"
            f"{response.text[:2000]}"
        )
    payload = response.json()
    if not isinstance(payload, dict):
        raise RuntimeError(f"Unexpected GET payload for {identifier}")
    return payload


def link_accessions(value: Any, prefix: str) -> set[str]:
    if value is None:
        items = []
    elif isinstance(value, (list, tuple, set)):
        items = list(value)
    else:
        items = [value]

    output: set[str] = set()
    for item in items:
        candidates: list[str] = []
        if isinstance(item, dict):
            for key in ("accession", "@id", "uuid"):
                if item.get(key):
                    candidates.append(str(item[key]))
            aliases = item.get("aliases", [])
            if isinstance(aliases, str):
                aliases = [aliases]
            candidates.extend(map(str, aliases))
        else:
            candidates.append(str(item))

        for candidate in candidates:
            match = re.search(rf"\b{re.escape(prefix)}[A-Z0-9]+\b", candidate)
            if match:
                output.add(match.group(0))
    return output


def patch_portal_record(
    session: requests.Session,
    identifier: str,
    payload: dict,
) -> dict:
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    response = session.patch(
        portal_record_url(identifier),
        json=payload,
        headers=headers,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )
    if not response.ok:
        raise RuntimeError(
            f"PATCH failed for {identifier}: HTTP {response.status_code}\n"
            f"{response.text[:4000]}"
        )
    value = response.json()
    if isinstance(value, dict) and isinstance(value.get("@graph"), list) and value["@graph"]:
        return value["@graph"][0]
    if isinstance(value, dict):
        return value
    raise RuntimeError(f"Unexpected PATCH response for {identifier}")
```

## Markdown cell 5

## 3. Materialize the embedded reference templates and submission map

## Code cell 6

```python
reference_paths = {}
for modality, text in CORRECTED_REFERENCE_TEXT.items():
    yaml_path = REFERENCE_DIR / f"corrected_{modality}_reference.yaml"
    gz_path = REFERENCE_DIR / f"corrected_{modality}_reference.yaml.gz"
    yaml_path.write_text(text, encoding="utf-8")
    write_deterministic_gzip(text.encode("utf-8"), gz_path)
    reference_paths[modality] = yaml_path

submission_map = pd.DataFrame(SUBMISSION_RECORDS)
submission_map["i"] = pd.to_numeric(submission_map["i"], errors="raise").astype(int)
submission_map = submission_map.sort_values("i").reset_index(drop=True)

if len(submission_map) != 48:
    raise RuntimeError(f"Expected 48 submission rows, observed {len(submission_map)}")
if submission_map["alias"].duplicated().any():
    raise RuntimeError("Duplicate ConfigurationFile aliases in embedded mapping")
if submission_map["accession"].duplicated().any():
    raise RuntimeError("Duplicate ConfigurationFile accessions in embedded mapping")
if submission_map.groupby("modality").size().to_dict() != {"GEX": 24, "gRNA": 24}:
    raise RuntimeError("Expected exactly 24 GEX and 24 gRNA records")

embedded_map_path = REPORT_DIR / "embedded_submission_mapping.csv"
submission_map.to_csv(embedded_map_path, index=False, lineterminator="\n")

print(submission_map.groupby("modality").size())
print("Reference templates:")
for modality, path in reference_paths.items():
    print(f"  {modality}: {path}")
print("Embedded map:", embedded_map_path)
```

### Text output (stream)

```text
modality
GEX     24
gRNA    24
dtype: int64
Reference templates:
  GEX: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reference_templates/corrected_GEX_reference.yaml
  gRNA: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reference_templates/corrected_gRNA_reference.yaml
Embedded map: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/embedded_submission_mapping.csv

```

## Markdown cell 7

## 4. Inventory and verify all 48 original submitted files

## Code cell 8

```python
source_rows: list[dict] = []
old_specs: dict[str, TaggedMap] = {}

for row in submission_map.to_dict("records"):
    basename = Path(row["submitted_file_name"]).name
    source_path = SOURCE_ROOT / basename

    if not source_path.is_file():
        raise FileNotFoundError(
            f"Missing original submitted file: {source_path}"
        )

    observed_md5 = md5_file(source_path)
    observed_size = source_path.stat().st_size

    if observed_md5 != row["md5sum"]:
        raise RuntimeError(
            f"Original source MD5 mismatch for {basename}: "
            f"expected {row['md5sum']}, observed {observed_md5}"
        )
    if int(observed_size) != int(row["file_size"]):
        raise RuntimeError(
            f"Original source size mismatch for {basename}: "
            f"expected {row['file_size']}, observed {observed_size}"
        )

    spec = load_tagged_yaml_gz(source_path)
    old_specs[row["accession"]] = spec

    expected_seqspec_of = {
        token.strip()
        for token in row["seqspec_of"].split(",")
        if token.strip()
    }
    observed_seqspec_of = set(all_file_ids(spec))
    if observed_seqspec_of != expected_seqspec_of:
        raise RuntimeError(
            f"Source YAML FASTQ IDs differ from seqspec_of for {row['accession']}: "
            f"YAML={sorted(observed_seqspec_of)}, mapping={sorted(expected_seqspec_of)}"
        )

    file_set_accession = file_set_accession_from_read_ids(spec)
    expected_file_set = re.search(r"IGVFDS[A-Z0-9]+", row["file_set"])
    if not expected_file_set or file_set_accession != expected_file_set.group(0):
        raise RuntimeError(
            f"Read IDs do not match file_set for {row['accession']}"
        )

    backup_path = BACKUP_DIR / basename
    if not backup_path.is_file() or md5_file(backup_path) != observed_md5:
        shutil.copy2(source_path, backup_path)

    source_rows.append(
        {
            "accession": row["accession"],
            "alias": row["alias"],
            "modality": row["modality"],
            "source_path": str(source_path),
            "source_md5": observed_md5,
            "source_size": observed_size,
            "original_date": str(spec["date"]),
            "file_set_accession": file_set_accession,
            "seqspec_of_count": len(observed_seqspec_of),
            "status": "PASS",
        }
    )

source_inventory = pd.DataFrame(source_rows)
source_inventory_path = REPORT_DIR / "original_source_inventory.tsv"
write_tsv(source_inventory, source_inventory_path)

print(source_inventory.groupby(["modality", "status"]).size())
print("Original source inventory:", source_inventory_path)
```

### Text output (stream)

```text
modality  status
GEX       PASS      24
gRNA      PASS      24
dtype: int64
Original source inventory: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/original_source_inventory.tsv

```

## Markdown cell 9

## 5. Generate the 48 fully corrected YAML and YAML.GZ files

## Code cell 10

```python
corrected_templates = {
    modality: load_tagged_yaml_text(text)
    for modality, text in CORRECTED_REFERENCE_TEXT.items()
}
corrected_template_fingerprints = {
    modality: normalized_template_object(template)
    for modality, template in corrected_templates.items()
}

generation_rows: list[dict] = []
diff_rows: list[dict] = []
corrected_specs: dict[str, TaggedMap] = {}

for row in submission_map.to_dict("records"):
    accession = row["accession"]
    modality = row["modality"]
    old_spec = old_specs[accession]
    corrected = build_corrected_spec(
        old_spec,
        corrected_templates[modality],
    )

    # Exact authoritative structural-template match after normalizing only the
    # explicitly preserved fields.
    observed_fingerprint = normalized_template_object(corrected)
    if observed_fingerprint != corrected_template_fingerprints[modality]:
        raise RuntimeError(
            f"Corrected {modality} structure differs from authoritative template: "
            f"{accession}"
        )

    # Explicit preservation checks.
    if corrected["date"] != old_spec["date"]:
        raise RuntimeError(f"Original date was not preserved for {accession}")
    if corrected["assay_id"] != old_spec["assay_id"]:
        raise RuntimeError(f"Original assay_id was not preserved for {accession}")

    old_reads = {read_name(read): read for read in old_spec["sequence_spec"]}
    new_reads = {read_name(read): read for read in corrected["sequence_spec"]}
    for key in ("read 1", "read 2"):
        if new_reads[key]["read_id"] != old_reads[key]["read_id"]:
            raise RuntimeError(f"read_id changed for {accession} {key}")
        if new_reads[key]["files"] != old_reads[key]["files"]:
            raise RuntimeError(f"FASTQ file metadata changed for {accession} {key}")

    expected_seqspec_of = {
        token.strip()
        for token in row["seqspec_of"].split(",")
        if token.strip()
    }
    if set(all_file_ids(corrected)) != expected_seqspec_of:
        raise RuntimeError(f"Corrected seqspec_of mismatch for {accession}")

    assert_joined_region_consistency(corrected)

    basename_gz = Path(row["submitted_file_name"]).name
    if not basename_gz.endswith(".yaml.gz"):
        raise RuntimeError(f"Unexpected source filename: {basename_gz}")
    basename_yaml = basename_gz[:-3]

    yaml_path = YAML_DIR / basename_yaml
    gz_path = YAML_GZ_DIR / basename_gz

    yaml_text = dump_tagged_yaml(corrected)
    yaml_bytes = yaml_text.encode("utf-8")
    yaml_path.write_bytes(yaml_bytes)
    write_deterministic_gzip(yaml_bytes, gz_path)

    # Disk round trip through both plain YAML and compressed YAML.
    if load_tagged_yaml_text(yaml_path.read_text(encoding="utf-8")) != corrected:
        raise RuntimeError(f"Plain YAML round trip failed for {accession}")
    if load_tagged_yaml_gz(gz_path) != corrected:
        raise RuntimeError(f"Gzip YAML round trip failed for {accession}")

    corrected_specs[accession] = corrected

    new_description_prefix = row["description"].split("; generated using", 1)[0]
    reference_id = "IGVFFI6906STSO" if modality == "GEX" else "IGVFFI4006RNQZ"
    new_description = (
        f"{new_description_prefix}; regenerated using the corrected {modality} "
        f"structural template represented by {reference_id}, while preserving "
        "the original sample-specific read/file metadata and date."
    )

    for change in recursive_diff(old_spec, corrected):
        diff_rows.append(
            {
                "accession": accession,
                "modality": modality,
                **change,
            }
        )

    generation_rows.append(
        {
            "accession": accession,
            "alias": row["alias"],
            "modality": modality,
            "file_set": row["file_set"],
            "seqspec_of": row["seqspec_of"],
            "original_date": str(old_spec["date"]),
            "source_md5": row["md5sum"],
            "source_size": int(row["file_size"]),
            "corrected_yaml_path": str(yaml_path),
            "corrected_gz_path": str(gz_path),
            "corrected_md5": md5_file(gz_path),
            "corrected_sha256": sha256_file(gz_path),
            "corrected_size": gz_path.stat().st_size,
            "new_description": new_description,
            "read2_pattern": expected_read2_pattern(corrected),
            "date_preserved": corrected["date"] == old_spec["date"],
            "read_metadata_preserved": True,
            "template_fingerprint_match": True,
        }
    )

generation_manifest = pd.DataFrame(generation_rows)
generation_manifest_path = REPORT_DIR / "corrected_generation_manifest.tsv"
write_tsv(generation_manifest, generation_manifest_path)

change_report = pd.DataFrame(diff_rows)
change_report_path = REPORT_DIR / "approved_structural_change_report.tsv"
write_tsv(change_report, change_report_path)

if generation_manifest["corrected_md5"].eq(generation_manifest["source_md5"]).any():
    raise RuntimeError("At least one corrected compressed file has the old MD5")

print(generation_manifest.groupby("modality").agg(files=("accession", "size"), total_bytes=("corrected_size", "sum")))
print("Generation manifest:", generation_manifest_path)
print("Change report:", change_report_path)
```

### Text output (stream)

```text
          files  total_bytes
modality                    
GEX          24        31995
gRNA         24        35136
Generation manifest: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/corrected_generation_manifest.tsv
Change report: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/approved_structural_change_report.tsv

```

## Markdown cell 11

## 6. Structural and Read 2 validation

## Code cell 12

```python
EXPECTED_READ2_PATTERN = (
    "XXXXXXXXXX"
    "NNNNNNNN"
    "GTGGCCGATGTTTCGCATCGGCGTACGACT"
    "NNNNNNNN"
    "ATCCACGTGCTTGAGACTGTGG"
    "NNNNNNNN"
)

EXPECTED_GEX_REGION_ORDER = [
    "illumina_p5",
    "truseq_read1",
    "cdna",
    "poly_A",
    "barcode-1",
    "linker-1",
    "barcode-2",
    "linker-2",
    "barcode-3",
    "umi",
    "truseq_read2",
    "illumina_p7",
]

EXPECTED_GRNA_REGION_ORDER = [
    "illumina_p5",
    "truseq_read1",
    "random_sequence",
    "prefix",
    "sgrna_target",
    "suffix",
    "barcode-1",
    "linker-1",
    "barcode-2",
    "linker-2",
    "barcode-3",
    "umi",
    "truseq_read2",
    "illumina_p7",
]

validation_rows: list[dict] = []

for row in generation_manifest.to_dict("records"):
    spec = corrected_specs[row["accession"]]
    modality = row["modality"]
    roots = spec["library_spec"]
    root = roots[0]
    order = [str(region["region_id"]) for region in root["regions"]]
    regions = region_by_id(spec)

    read2 = {
        read_name(read): read
        for read in spec["sequence_spec"]
    }["read 2"]

    checks = {
        "read2_negative_strand": read2["strand"] == "neg",
        "read2_primer": read2["primer_id"] == "truseq_read2",
        "read2_length_86": int(read2["min_len"]) == 86 and int(read2["max_len"]) == 86,
        "read2_pattern": expected_read2_pattern(spec) == EXPECTED_READ2_PATTERN,
        "linker1_physical_sequence": regions["linker-1"]["sequence"] == "CCACAGTCTCAAGCACGTGGAT",
        "linker2_physical_sequence": regions["linker-2"]["sequence"] == "AGTCGTACGCCGATGCGAAACATCGGCCAC",
        "barcode2_md5": regions["barcode-2"]["onlist"]["md5"] == "7d310eca0f1030f13fec3395d044142d",
        "barcode3_md5": regions["barcode-3"]["onlist"]["md5"] == "7d310eca0f1030f13fec3395d044142d",
        "barcode_onlists_remote": all(
            regions[key]["onlist"].get("location") == "remote"
            for key in ("barcode-1", "barcode-2", "barcode-3")
        ),
    }

    if modality == "GEX":
        checks.update(
            {
                "region_order": order == EXPECTED_GEX_REGION_ORDER,
                "poly_A": regions["poly_A"]["sequence"] == "AAAAAAAAAAAAAA",
                "corrected_sequence_protocol": spec["sequence_protocol"] == "Illumina NovaSeq X (NTR:0000765)",
            }
        )
    else:
        checks.update(
            {
                "region_order": order == EXPECTED_GRNA_REGION_ORDER,
                "corrected_read1_primer": regions["truseq_read1"]["sequence"] == "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                "corrected_read2_primer": regions["truseq_read2"]["sequence"] == "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                "sgrna_target_semantics": (
                    regions["sgrna_target"]["region_type"] == "gdna"
                    and regions["sgrna_target"]["name"] == "sgrna_target"
                ),
                "guide_onlist_remote": regions["sgrna_target"]["onlist"].get("location") == "remote",
                "corrected_sequence_protocol": spec["sequence_protocol"] == "Illumina NovaSeq X (NTR:0000765)",
            }
        )

    failed_checks = sorted(key for key, value in checks.items() if not value)
    validation_rows.append(
        {
            "accession": row["accession"],
            "modality": modality,
            "checks": len(checks),
            "failed_checks": ",".join(failed_checks),
            "status": "PASS" if not failed_checks else "FAIL",
        }
    )

structural_validation = pd.DataFrame(validation_rows)
structural_validation_path = REPORT_DIR / "structural_validation.tsv"
write_tsv(structural_validation, structural_validation_path)

print(structural_validation.groupby(["modality", "status"]).size())
if not structural_validation["status"].eq("PASS").all():
    raise RuntimeError(
        "At least one corrected seqspec failed structural validation. "
        f"Review {structural_validation_path}"
    )
print("Structural validation:", structural_validation_path)
```

### Text output (stream)

```text
modality  status
GEX       PASS      24
gRNA      PASS      24
dtype: int64
Structural validation: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/structural_validation.tsv

```

## Markdown cell 13

## 7. Validate all 48 corrected YAMLs with seqspec

## Code cell 14

```python
# Replacement Section 7: authenticated validation of private IGVF FASTQ files.
#
# This cell assumes Sections 1-6 of resubmit_corrected_300gene_seqspecs.ipynb
# have already been run, producing generation_manifest and corrected YAML files.

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


def seqspec_subprocess_env() -> dict[str, str]:
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    mpl_dir = WORK_DIR / "matplotlib_config"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    environment["MPLCONFIGDIR"] = str(mpl_dir)
    return environment


def find_seqspec_executable() -> Path:
    path_candidate = shutil.which("seqspec")
    candidates = [
        *SEQSPEC_EXECUTABLE_CANDIDATES,
        Path(path_candidate) if path_candidate else None,
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "No seqspec executable was found. Expected one of: "
        + ", ".join(map(str, SEQSPEC_EXECUTABLE_CANDIDATES))
    )


def record_type_names(record: dict) -> set[str]:
    value = record.get("@type", [])
    if isinstance(value, str):
        return {value}
    if isinstance(value, list):
        return {str(item) for item in value}
    return set()


def expected_sequence_file_download_url(
    file_id: str,
    filename: str,
) -> str:
    return (
        API_BASE.rstrip("/")
        + "/sequence-files/"
        + file_id
        + "/@@download/"
        + filename
    )


def get_portal_record_with_retry(
    session,
    identifier: str,
    *,
    attempts: int = 3,
) -> dict:
    last_error: Exception | None = None

    for attempt in range(1, attempts + 1):
        try:
            return get_portal_record(
                session,
                identifier,
                frame="object",
                database=True,
            )
        except Exception as error:
            last_error = error
            if attempt < attempts:
                time.sleep(min(2 ** (attempt - 1), 5))

    assert last_error is not None
    raise last_error


def authenticated_range_probe(
    session,
    url: str,
    *,
    attempts: int = 3,
) -> dict:
    """Verify a private download URL without downloading the FASTQ body."""
    last_error = ""
    last_status: int | str = ""
    final_host = ""

    for attempt in range(1, attempts + 1):
        try:
            response = session.get(
                url,
                headers={
                    "Range": "bytes=0-0",
                    "Accept": "*/*",
                },
                stream=True,
                allow_redirects=True,
                timeout=REQUEST_TIMEOUT_SECONDS,
            )
            last_status = int(response.status_code)
            final_host = urlparse(str(response.url)).hostname or ""
            success = response.status_code in {200, 206}
            response.close()

            if success:
                return {
                    "download_probe_success": True,
                    "download_probe_status": last_status,
                    "download_final_host": final_host,
                    "download_probe_error": "",
                }

            last_error = f"HTTP {last_status}"
        except Exception as error:
            last_error = f"{type(error).__name__}: {error}"

        if attempt < attempts:
            time.sleep(min(2 ** (attempt - 1), 5))

    return {
        "download_probe_success": False,
        "download_probe_status": last_status,
        "download_final_host": final_host,
        "download_probe_error": last_error,
    }


def extract_seqspec_fastq_objects(spec: dict) -> list[dict]:
    rows: list[dict] = []

    for read in spec.get("sequence_spec", []):
        read_id = str(read.get("read_id", ""))
        read_label = str(read.get("name", ""))

        for file_object in read.get("files", []) or []:
            file_id = str(file_object.get("file_id", "")).strip()
            filename = str(file_object.get("filename", "")).strip()
            url = str(file_object.get("url", "")).strip()
            urltype = str(file_object.get("urltype", "")).strip()
            filetype = str(file_object.get("filetype", "")).strip()
            md5 = str(file_object.get("md5", "")).strip()

            try:
                filesize = int(file_object.get("filesize"))
            except (TypeError, ValueError):
                filesize = -1

            if not re.fullmatch(r"IGVFFI[A-Z0-9]+", file_id):
                raise RuntimeError(
                    f"Invalid SequenceFile accession in {read_id}: {file_id!r}"
                )

            rows.append(
                {
                    "read_id": read_id,
                    "read_label": read_label,
                    "file_id": file_id,
                    "filename": filename,
                    "url": url,
                    "urltype": urltype,
                    "filetype": filetype,
                    "filesize": filesize,
                    "md5": md5,
                }
            )

    return rows


def classify_seqspec_error_lines(
    error_lines: list[str],
    expected_filenames: set[str],
) -> tuple[list[str], list[str]]:
    """Separate unauthenticated private-file errors from real seqspec errors."""
    expected_remote_errors: list[str] = []
    unexpected_errors: list[str] = []

    pattern = re.compile(
        r"^\[error(?:\s+\d+)?\]\s+(.+?)\s+does not exist$",
        flags=re.IGNORECASE,
    )

    for line in error_lines:
        match = pattern.fullmatch(line.strip())
        if match and match.group(1) in expected_filenames:
            expected_remote_errors.append(line.strip())
        else:
            unexpected_errors.append(line.strip())

    return expected_remote_errors, unexpected_errors


# ---------------------------------------------------------------------------
# 1. Locate the seqspec environment and prepare one batch validator process.
# ---------------------------------------------------------------------------
SEQSPEC_EXECUTABLE = find_seqspec_executable()
SEQSPEC_PYTHON = SEQSPEC_EXECUTABLE.parent / "python"

if not SEQSPEC_PYTHON.is_file():
    try:
        first_line = SEQSPEC_EXECUTABLE.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()[0]
    except Exception:
        first_line = ""

    if first_line.startswith("#!"):
        shebang_python = Path(first_line[2:].strip())
        if shebang_python.is_file():
            SEQSPEC_PYTHON = shebang_python

if not SEQSPEC_PYTHON.is_file():
    raise FileNotFoundError(
        "Could not identify the Python executable for the selected seqspec "
        f"command: {SEQSPEC_EXECUTABLE}"
    )


# ---------------------------------------------------------------------------
# 2. Load all FASTQ references from the corrected YAMLs.
# ---------------------------------------------------------------------------
reference_rows: list[dict] = []
expected_filenames_by_seqspec: dict[str, set[str]] = {}
file_metadata_by_accession: dict[str, dict] = {}
file_seqspec_membership: dict[str, set[str]] = {}

for manifest_row in generation_manifest.to_dict("records"):
    config_accession = str(manifest_row["accession"])
    yaml_path = Path(manifest_row["corrected_yaml_path"])
    spec = load_tagged_yaml_text(
        yaml_path.read_text(encoding="utf-8")
    )
    file_objects = extract_seqspec_fastq_objects(spec)

    if len(file_objects) != 2:
        raise RuntimeError(
            f"{config_accession} must reference exactly two FASTQs; "
            f"observed {len(file_objects)}"
        )

    expected_filenames_by_seqspec[config_accession] = {
        item["filename"] for item in file_objects
    }

    for item in file_objects:
        file_id = item["file_id"]
        candidate_metadata = {
            key: item[key]
            for key in (
                "filename",
                "url",
                "urltype",
                "filetype",
                "filesize",
                "md5",
            )
        }

        if (
            file_id in file_metadata_by_accession
            and file_metadata_by_accession[file_id] != candidate_metadata
        ):
            raise RuntimeError(
                f"SequenceFile {file_id} has inconsistent metadata across YAMLs"
            )

        file_metadata_by_accession[file_id] = candidate_metadata
        file_seqspec_membership.setdefault(file_id, set()).add(config_accession)

        reference_rows.append(
            {
                "configuration_file_accession": config_accession,
                "modality": manifest_row["modality"],
                "yaml_path": str(yaml_path),
                **item,
            }
        )

if len(file_metadata_by_accession) != 96:
    raise RuntimeError(
        "Expected 96 unique SequenceFiles across 48 seqspecs; "
        f"observed {len(file_metadata_by_accession)}"
    )

seqspec_fastq_references = pd.DataFrame(reference_rows)
seqspec_fastq_references_path = (
    REPORT_DIR / "seqspec_fastq_references.tsv"
)
write_tsv(
    seqspec_fastq_references,
    seqspec_fastq_references_path,
)


# ---------------------------------------------------------------------------
# 3. Authenticate and verify all 96 private SequenceFiles.
# ---------------------------------------------------------------------------
enter_igvf_credentials()
portal_session = make_portal_session()

portal_file_rows: list[dict] = []
portal_file_pass_by_accession: dict[str, bool] = {}

for index, file_id in enumerate(
    sorted(file_metadata_by_accession),
    start=1,
):
    expected = file_metadata_by_accession[file_id]
    record = get_portal_record_with_retry(
        portal_session,
        file_id,
    )

    types = record_type_names(record)
    observed_accession = str(record.get("accession", "")).strip()
    observed_md5 = str(record.get("md5sum", "")).strip()
    observed_file_format = str(record.get("file_format", "")).strip()

    try:
        observed_size = int(record.get("file_size"))
    except (TypeError, ValueError):
        observed_size = -1

    expected_filename = f"{file_id}.fastq.gz"
    expected_url = expected_sequence_file_download_url(
        file_id,
        expected_filename,
    )

    checks = {
        "record_type": "SequenceFile" in types,
        "accession": observed_accession == file_id,
        "file_format": observed_file_format == "fastq",
        "yaml_filename": expected["filename"] == expected_filename,
        "yaml_filetype": expected["filetype"] == "fastq",
        "yaml_urltype": expected["urltype"] == "https",
        "yaml_url": expected["url"] == expected_url,
        "md5": observed_md5 == expected["md5"],
        "file_size": observed_size == int(expected["filesize"]),
    }

    probe = authenticated_range_probe(
        portal_session,
        expected["url"],
    )
    checks["authenticated_download_probe"] = bool(
        probe["download_probe_success"]
    )

    failed_checks = sorted(
        key for key, passed in checks.items() if not passed
    )
    passed = not failed_checks
    portal_file_pass_by_accession[file_id] = passed

    portal_file_rows.append(
        {
            "sequence_file_accession": file_id,
            "configuration_file_count": len(
                file_seqspec_membership[file_id]
            ),
            "configuration_file_accessions": ",".join(
                sorted(file_seqspec_membership[file_id])
            ),
            "portal_status": str(record.get("status", "")),
            "portal_upload_status": str(record.get("upload_status", "")),
            "expected_filename": expected["filename"],
            "expected_md5": expected["md5"],
            "portal_md5": observed_md5,
            "expected_file_size": int(expected["filesize"]),
            "portal_file_size": observed_size,
            **probe,
            "failed_checks": ",".join(failed_checks),
            "status": "PASS" if passed else "FAIL",
        }
    )

    if index == 1 or index % 12 == 0 or index == len(file_metadata_by_accession):
        print(
            f"Authenticated SequenceFile validation "
            f"{index}/{len(file_metadata_by_accession)}",
            flush=True,
        )

portal_sequence_file_validation = pd.DataFrame(portal_file_rows)
portal_sequence_file_validation_path = (
    REPORT_DIR / "seqspec_sequence_file_portal_validation.tsv"
)
write_tsv(
    portal_sequence_file_validation,
    portal_sequence_file_validation_path,
)

print("\nAuthenticated SequenceFile validation:")
print(
    portal_sequence_file_validation["status"].value_counts(
        dropna=False
    )
)

if not portal_sequence_file_validation["status"].eq("PASS").all():
    raise RuntimeError(
        "At least one referenced SequenceFile failed authenticated Portal "
        f"validation. Review {portal_sequence_file_validation_path}"
    )


# ---------------------------------------------------------------------------
# 4. Run all 48 seqspec checks in one external Python process.
# ---------------------------------------------------------------------------
batch_input_path = WORK_DIR / "seqspec_batch_validation_input.json"
batch_output_path = WORK_DIR / "seqspec_batch_validation_output.json"
batch_runner_path = WORK_DIR / "seqspec_batch_validator.py"
batch_log_path = LOG_DIR / "seqspec_batch_validation.log"

batch_input = [
    {
        "accession": str(row["accession"]),
        "modality": str(row["modality"]),
        "yaml_path": str(row["corrected_yaml_path"]),
    }
    for row in generation_manifest.to_dict("records")
]
batch_input_path.write_text(
    json.dumps(batch_input, indent=2) + "\n",
    encoding="utf-8",
)

batch_runner_source = r'''from __future__ import annotations

import json
import sys
from pathlib import Path

from seqspec.seqspec_check import seqspec_check
from seqspec.utils import load_spec


def main() -> None:
    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    rows = json.loads(input_path.read_text(encoding="utf-8"))
    results = []

    for row in rows:
        try:
            spec = load_spec(Path(row["yaml_path"]), strict=False)
            errors = seqspec_check(spec, "igvf_onlist_skip")
            results.append(
                {
                    **row,
                    "returncode": 0,
                    "exception": "",
                    "errors": errors,
                }
            )
        except Exception as error:
            results.append(
                {
                    **row,
                    "returncode": 1,
                    "exception": f"{type(error).__name__}: {error}",
                    "errors": [],
                }
            )

    output_path.write_text(
        json.dumps(results, indent=2, default=str) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
'''
batch_runner_path.write_text(
    batch_runner_source,
    encoding="utf-8",
)

compile_result = subprocess.run(
    [
        str(SEQSPEC_PYTHON),
        "-m",
        "py_compile",
        str(batch_runner_path),
    ],
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    env=seqspec_subprocess_env(),
    check=False,
)
if compile_result.returncode != 0:
    raise RuntimeError(
        "The seqspec batch validator failed syntax validation:\n"
        + compile_result.stdout
    )

batch_result = subprocess.run(
    [
        str(SEQSPEC_PYTHON),
        str(batch_runner_path),
        str(batch_input_path),
        str(batch_output_path),
    ],
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    env=seqspec_subprocess_env(),
    check=False,
)
batch_log_path.write_text(
    batch_result.stdout,
    encoding="utf-8",
)

if batch_result.returncode != 0 or not batch_output_path.is_file():
    raise RuntimeError(
        "Seqspec batch validation process failed. Review "
        f"{batch_log_path}"
    )

batch_results = json.loads(
    batch_output_path.read_text(encoding="utf-8")
)
if not isinstance(batch_results, list) or len(batch_results) != 48:
    raise RuntimeError(
        "Seqspec batch validator did not return exactly 48 results"
    )

batch_by_accession = {
    str(row["accession"]): row
    for row in batch_results
}
if set(batch_by_accession) != set(generation_manifest["accession"].astype(str)):
    raise RuntimeError(
        "Seqspec batch validator result accessions do not match the manifest"
    )


# ---------------------------------------------------------------------------
# 5. Accept only the known unauthenticated private-file errors after the
#    authenticated record and download verification has passed.
# ---------------------------------------------------------------------------
validation_rows: list[dict] = []

for index, row in enumerate(
    generation_manifest.to_dict("records"),
    start=1,
):
    config_accession = str(row["accession"])
    yaml_path = Path(row["corrected_yaml_path"])
    result = batch_by_accession[config_accession]
    error_objects = result.get("errors", [])

    expected_remote_errors: list[str] = []
    unexpected_errors: list[str] = []

    for error in error_objects:
        error_type = str(error.get("error_type", ""))
        error_object = str(error.get("error_object", ""))
        error_message = str(error.get("error_message", ""))

        match = re.fullmatch(
            r"(.+?)\s+does not exist",
            error_message,
            flags=re.IGNORECASE,
        )

        if (
            error_type == "check_read_files_exist"
            and error_object == "file"
            and match
            and match.group(1) in expected_filenames_by_seqspec[config_accession]
        ):
            expected_remote_errors.append(error_message)
        else:
            unexpected_errors.append(
                json.dumps(error, sort_keys=True, default=str)
            )

    referenced_file_ids = sorted(
        seqspec_fastq_references.loc[
            seqspec_fastq_references[
                "configuration_file_accession"
            ].eq(config_accession),
            "file_id",
        ].unique()
    )
    portal_files_passed = all(
        portal_file_pass_by_accession[file_id]
        for file_id in referenced_file_ids
    )

    exception = str(result.get("exception", ""))
    passed = (
        int(result.get("returncode", 1)) == 0
        and not exception
        and not unexpected_errors
        and portal_files_passed
    )

    validation_rows.append(
        {
            "accession": config_accession,
            "modality": row["modality"],
            "yaml_path": str(yaml_path),
            "returncode": int(result.get("returncode", 1)),
            "exception": exception,
            "referenced_sequence_files": ",".join(referenced_file_ids),
            "portal_sequence_files_passed": portal_files_passed,
            "expected_private_url_errors": " | ".join(
                expected_remote_errors
            ),
            "unexpected_errors": " | ".join(unexpected_errors),
            "validation_basis": (
                "seqspec Python API structural checks + authenticated IGVF "
                "record, metadata, canonical URL, and ranged-download verification"
            ),
            "batch_log_path": str(batch_log_path),
            "status": "PASS" if passed else "FAIL",
        }
    )

    if index == 1 or index % 12 == 0 or index == len(generation_manifest):
        print(
            f"Validated corrected seqspec {index}/{len(generation_manifest)}",
            flush=True,
        )

seqspec_validation = pd.DataFrame(validation_rows)
seqspec_validation_path = REPORT_DIR / "seqspec_cli_validation.tsv"
write_tsv(seqspec_validation, seqspec_validation_path)

print("\nCorrected seqspec validation:")
print(seqspec_validation.groupby(["modality", "status"]).size())

if not seqspec_validation["status"].eq("PASS").all():
    raise RuntimeError(
        "At least one corrected seqspec failed validation. Review "
        f"{seqspec_validation_path}"
    )

print("\nSeqspec validation:", seqspec_validation_path)
print(
    "Authenticated SequenceFile validation:",
    portal_sequence_file_validation_path,
)

```

### Text output (stream)

```text
IGVF credentials are set for the current kernel.
Authenticated SequenceFile validation 1/96
Authenticated SequenceFile validation 12/96
Authenticated SequenceFile validation 24/96
Authenticated SequenceFile validation 36/96
Authenticated SequenceFile validation 48/96
Authenticated SequenceFile validation 60/96
Authenticated SequenceFile validation 72/96
Authenticated SequenceFile validation 84/96
Authenticated SequenceFile validation 96/96

Authenticated SequenceFile validation:
status
PASS    96
Name: count, dtype: int64
Validated corrected seqspec 1/48
Validated corrected seqspec 12/48
Validated corrected seqspec 24/48
Validated corrected seqspec 36/48
Validated corrected seqspec 48/48

Corrected seqspec validation:
modality  status
GEX       PASS      24
gRNA      PASS      24
dtype: int64

Seqspec validation: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/seqspec_cli_validation.tsv
Authenticated SequenceFile validation: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/seqspec_sequence_file_portal_validation.tsv

```

## Code cell 15

```python
def seqspec_subprocess_env() -> dict[str, str]:
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    mpl_dir = WORK_DIR / "matplotlib_config"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    environment["MPLCONFIGDIR"] = str(mpl_dir)
    return environment


def find_seqspec_executable() -> Path:
    path_candidate = shutil.which("seqspec")
    candidates = [
        *SEQSPEC_EXECUTABLE_CANDIDATES,
        Path(path_candidate) if path_candidate else None,
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "No seqspec executable was found. Expected one of: "
        + ", ".join(map(str, SEQSPEC_EXECUTABLE_CANDIDATES))
    )


SEQSPEC_EXECUTABLE = find_seqspec_executable()
help_result = subprocess.run(
    [str(SEQSPEC_EXECUTABLE), "check", "--help"],
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    env=seqspec_subprocess_env(),
    check=False,
)
if help_result.returncode != 0:
    raise RuntimeError(
        f"seqspec check --help failed:\n{help_result.stdout}"
    )

supports_skip = "--skip" in help_result.stdout
validation_rows = []

for index, row in enumerate(generation_manifest.to_dict("records"), start=1):
    yaml_path = Path(row["corrected_yaml_path"])
    log_path = LOG_DIR / f"{yaml_path.stem}_seqspec_check.log"

    command = [str(SEQSPEC_EXECUTABLE), "check"]
    if supports_skip:
        command.extend(["--skip", "igvf_onlist_skip"])
    command.append(str(yaml_path))

    result = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=seqspec_subprocess_env(),
        check=False,
    )
    log_path.write_text(result.stdout, encoding="utf-8")

    error_lines = [
        line.strip()
        for line in result.stdout.splitlines()
        if re.search(r"\[error(?:\s+\d+)?\]", line, flags=re.IGNORECASE)
    ]
    invocation_error = (
        result.returncode == 2
        or re.search(
            r"usage:|unrecognized arguments|invalid choice",
            result.stdout,
            flags=re.IGNORECASE,
        )
        is not None
    )
    passed = result.returncode == 0 and not error_lines and not invocation_error

    validation_rows.append(
        {
            "accession": row["accession"],
            "modality": row["modality"],
            "yaml_path": str(yaml_path),
            "returncode": result.returncode,
            "error_lines": " | ".join(error_lines),
            "log_path": str(log_path),
            "status": "PASS" if passed else "FAIL",
        }
    )

    if index == 1 or index % 12 == 0 or index == len(generation_manifest):
        print(f"Validated {index}/{len(generation_manifest)}", flush=True)

seqspec_validation = pd.DataFrame(validation_rows)
seqspec_validation_path = REPORT_DIR / "seqspec_cli_validation.tsv"
write_tsv(seqspec_validation, seqspec_validation_path)

print(seqspec_validation.groupby(["modality", "status"]).size())
if not seqspec_validation["status"].eq("PASS").all():
    raise RuntimeError(
        "At least one corrected seqspec failed CLI validation. "
        f"Review {seqspec_validation_path}"
    )
print("seqspec validation:", seqspec_validation_path)
```

### Text output (stream)

```text
Validated 1/48
Validated 12/48
Validated 24/48
Validated 36/48
Validated 48/48
modality  status
GEX       PASS      24
gRNA      PASS      24
dtype: int64
seqspec validation: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/seqspec_cli_validation.tsv

```

## Markdown cell 16

## 8. Build a corrected resubmission package

## Code cell 17

```python
package_path = OUTPUT_ROOT / "corrected_300gene_seqspec_resubmission_package.zip"
with zipfile.ZipFile(
    package_path,
    mode="w",
    compression=zipfile.ZIP_DEFLATED,
    compresslevel=9,
) as archive:
    for directory in (YAML_DIR, YAML_GZ_DIR, REFERENCE_DIR, REPORT_DIR):
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(OUTPUT_ROOT))

package_manifest = {
    "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    "files": len(generation_manifest),
    "gex_files": int((generation_manifest["modality"] == "GEX").sum()),
    "grna_files": int((generation_manifest["modality"] == "gRNA").sum()),
    "package": str(package_path),
    "package_sha256": sha256_file(package_path),
}
package_manifest_path = REPORT_DIR / "package_manifest.json"
package_manifest_path.write_text(
    json.dumps(package_manifest, indent=2) + "\n",
    encoding="utf-8",
)

print(json.dumps(package_manifest, indent=2))
```

### Text output (stream)

```text
{
  "created_at": "2026-10-15T07:13:59-0700",
  "files": 48,
  "gex_files": 24,
  "grna_files": 24,
  "package": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/corrected_300gene_seqspec_resubmission_package.zip",
  "package_sha256": "d0de561dbf0ccc892e12bb794b923eec2d10b835f4be820123216f1fb7d4067b"
}

```

## Markdown cell 18

## 9. Enter IGVF credentials and run the Portal preflight

## Code cell 19

```python
enter_igvf_credentials()
portal_session = make_portal_session()

portal_rows: list[dict] = []
portal_records: dict[str, dict] = {}

blocked_statuses = {
    "released",
    "archived",
    "revoked",
    "replaced",
    "deleted",
}

manifest_by_accession = generation_manifest.set_index("accession").to_dict("index")
map_by_accession = submission_map.set_index("accession").to_dict("index")

for index, accession in enumerate(generation_manifest["accession"], start=1):
    record = get_portal_record(
        portal_session,
        accession,
        frame="object",
        database=True,
    )
    portal_records[accession] = record

    mapping = map_by_accession[accession]
    generated = manifest_by_accession[accession]

    aliases = record.get("aliases", [])
    if isinstance(aliases, str):
        aliases = [aliases]
    aliases = set(map(str, aliases))

    observed_file_set = link_accessions(record.get("file_set", ""), "IGVFDS")
    expected_file_set_match = re.search(r"IGVFDS[A-Z0-9]+", mapping["file_set"])
    expected_file_set = {expected_file_set_match.group(0)} if expected_file_set_match else set()

    observed_seqspec_of = link_accessions(record.get("seqspec_of", []), "IGVFFI")
    expected_seqspec_of = {
        value.strip()
        for value in mapping["seqspec_of"].split(",")
        if value.strip()
    }

    current_md5 = str(record.get("md5sum", ""))
    new_md5 = generated["corrected_md5"]

    if current_md5 == new_md5:
        action = "SKIP_ALREADY_CORRECT"
    elif current_md5 == mapping["md5sum"]:
        action = "UPDATE_REQUIRED"
    else:
        action = "BLOCK_UNEXPECTED_PORTAL_MD5"

    status = str(record.get("status", ""))
    errors = []
    if mapping["alias"] not in aliases:
        errors.append("alias mismatch")
    if observed_file_set != expected_file_set:
        errors.append("file_set mismatch")
    if observed_seqspec_of != expected_seqspec_of:
        errors.append("seqspec_of mismatch")
    if str(record.get("content_type", "")) != "seqspec":
        errors.append("content_type mismatch")
    if str(record.get("file_format", "")) != "yaml":
        errors.append("file_format mismatch")
    if status in blocked_statuses:
        errors.append(f"blocked status: {status}")

    if errors:
        action = "BLOCK_RECORD_MISMATCH"

    portal_rows.append(
        {
            "accession": accession,
            "alias": mapping["alias"],
            "modality": mapping["modality"],
            "status": status,
            "upload_status": str(record.get("upload_status", "")),
            "current_md5": current_md5,
            "expected_old_md5": mapping["md5sum"],
            "corrected_md5": new_md5,
            "current_file_size": record.get("file_size", ""),
            "corrected_file_size": generated["corrected_size"],
            "current_submitted_file_name": str(record.get("submitted_file_name", "")),
            "corrected_submitted_file_name": generated["corrected_gz_path"],
            "file_set_match": observed_file_set == expected_file_set,
            "seqspec_of_match": observed_seqspec_of == expected_seqspec_of,
            "errors": "; ".join(errors),
            "action": action,
        }
    )

    backup_path = PORTAL_BACKUP_DIR / f"{accession}.json"
    backup_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    if index == 1 or index % 12 == 0 or index == len(generation_manifest):
        print(f"Portal preflight {index}/{len(generation_manifest)}", flush=True)

portal_plan = pd.DataFrame(portal_rows)
portal_plan_path = REPORT_DIR / "portal_update_plan.tsv"
write_tsv(portal_plan, portal_plan_path)

print(portal_plan["action"].value_counts(dropna=False))
print(portal_plan.groupby(["modality", "action"]).size())

blocked = portal_plan.loc[portal_plan["action"].str.startswith("BLOCK_")]
if not blocked.empty:
    raise RuntimeError(
        "Portal preflight contains blocking records. Review "
        f"{portal_plan_path}"
    )

print("Portal update plan:", portal_plan_path)
```

### Text output (stream)

```text
IGVF credentials are set for the current kernel.
Portal preflight 1/48
Portal preflight 12/48
Portal preflight 24/48
Portal preflight 36/48
Portal preflight 48/48
action
UPDATE_REQUIRED    48
Name: count, dtype: int64
modality  action         
GEX       UPDATE_REQUIRED    24
gRNA      UPDATE_REQUIRED    24
dtype: int64
Portal update plan: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/portal_update_plan.tsv

```

## Markdown cell 20

## 10. Build metadata PATCH payloads and perform required dry runs

## Code cell 21

```python
update_operations: list[dict] = []

for row in generation_manifest.to_dict("records"):
    plan_row = portal_plan.loc[portal_plan["accession"].eq(row["accession"])].iloc[0]

    payload = {
        "record_id": row["accession"],
        "md5sum": row["corrected_md5"],
        "file_size": int(row["corrected_size"]),
        "submitted_file_name": row["corrected_gz_path"],
        "description": row["new_description"],
    }

    update_operations.append(
        {
            "accession": row["accession"],
            "modality": row["modality"],
            "action": plan_row["action"],
            "corrected_gz_path": row["corrected_gz_path"],
            "corrected_md5": row["corrected_md5"],
            "corrected_size": int(row["corrected_size"]),
            "payload": payload,
        }
    )

operations_path = WORK_DIR / "configuration_file_update_operations.json"
operations_path.write_text(
    json.dumps(update_operations, indent=2, ensure_ascii=False) + "\n",
    encoding="utf-8",
)

runner_path = WORK_DIR / "dry_run_configuration_file_patches.py"
runner_path.write_text(
    r"""from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from igvf_utils.connection import Connection


def main() -> None:
    operations_path = Path(sys.argv[1])
    operations = json.loads(operations_path.read_text(encoding="utf-8"))

    connection = Connection("prod", True)
    failures = []

    for index, operation in enumerate(operations, start=1):
        if operation["action"] == "SKIP_ALREADY_CORRECT":
            print(
                f"[{index}/{len(operations)}] SKIP {operation['accession']} "
                "already has the corrected MD5"
            )
            continue

        payload = dict(operation["payload"])
        accession = payload.pop("record_id")
        payload[connection.IGVFID_KEY] = accession

        print(f"[{index}/{len(operations)}] DRY RUN PATCH {accession}")
        try:
            connection.patch(payload, extend_array_values=False)
        except Exception as error:
            failures.append(
                f"{accession}: {type(error).__name__}: {error}"
            )

    if failures:
        raise RuntimeError("Dry-run failures:\n" + "\n".join(failures))


if __name__ == "__main__":
    main()
""",
    encoding="utf-8",
)

subprocess.run(
    [str(IGVF_PYTHON), "-m", "py_compile", str(runner_path)],
    check=True,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
)

DRY_RUN_MARKER = WORK_DIR / "METADATA_DRY_RUN_PASSED.json"

if not RUN_METADATA_DRY_RUNS:
    print(
        "Set RUN_METADATA_DRY_RUNS=True in Section 1 and rerun this cell. "
        "Production is blocked until the marker is created."
    )
else:
    environment = os.environ.copy()
    environment["IGVF_API_KEY"] = os.environ["IGVF_API_KEY"]
    environment["IGVF_SECRET_KEY"] = os.environ["IGVF_SECRET_KEY"]

    log_path = LOG_DIR / "configuration_file_metadata_patch_dryrun.log"
    result = subprocess.run(
        [str(IGVF_PYTHON), str(runner_path), str(operations_path)],
        cwd=str(LOG_DIR),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
        check=False,
    )
    log_path.write_text(result.stdout, encoding="utf-8")
    print(result.stdout)

    if result.returncode != 0 or re.search(
        r"Traceback|ValidationFailure|Failed to PATCH|Client Error",
        result.stdout,
        flags=re.IGNORECASE,
    ):
        raise RuntimeError(f"Metadata dry run failed. Review {log_path}")

    marker = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "generation_manifest_sha256": sha256_file(generation_manifest_path),
        "portal_plan_sha256": sha256_file(portal_plan_path),
        "operations_sha256": sha256_file(operations_path),
        "files": {
            row["accession"]: row["corrected_sha256"]
            for row in generation_manifest.to_dict("records")
        },
        "dry_run_log": str(log_path),
    }
    DRY_RUN_MARKER.write_text(
        json.dumps(marker, indent=2) + "\n",
        encoding="utf-8",
    )
    print("All metadata dry runs passed.")
    print("Dry-run marker:", DRY_RUN_MARKER)
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
2026-10-15 07:18:27,895:iu_debug:	DRY RUN is enabled.
2026-10-15 07:18:27,896:iu_debug:	submission=False: In non-submission mode.
[1/48] DRY RUN PATCH IGVFFI9366VIOI
2026-10-15 07:18:27,896:iu_debug:	
IN patch()
2026-10-15 07:18:28,097:iu_debug:	>>>>>>GET IGVFFI9366VIOI From DACC with URL https://api.data.igvf.org/IGVFFI9366VIOI/?format=json&frame=edit
2026-10-15 07:18:28,340:iu_debug:	<<<<<< PATCHING IGVFFI9366VIOI To IGVF database with URL https://api.data.igvf.org/IGVFFI9366VIOI and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D0_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "d8c1e16898161eec5bc237d77826d0fe",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:28,340:iu_debug:	DRY RUN is enabled.
[2/48] DRY RUN PATCH IGVFFI5296LCGA
2026-10-15 07:18:28,340:iu_debug:	
IN patch()
2026-10-15 07:18:28,340:iu_debug:	>>>>>>GET IGVFFI5296LCGA From DACC with URL https://api.data.igvf.org/IGVFFI5296LCGA/?format=json&frame=edit
2026-10-15 07:18:28,500:iu_debug:	<<<<<< PATCHING IGVFFI5296LCGA To IGVF database with URL https://api.data.igvf.org/IGVFFI5296LCGA and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D0_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "dcd154e3cdd84683e4dfab8e1a589962",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:28,501:iu_debug:	DRY RUN is enabled.
[3/48] DRY RUN PATCH IGVFFI3664MJLW
2026-10-15 07:18:28,501:iu_debug:	
IN patch()
2026-10-15 07:18:28,501:iu_debug:	>>>>>>GET IGVFFI3664MJLW From DACC with URL https://api.data.igvf.org/IGVFFI3664MJLW/?format=json&frame=edit
2026-10-15 07:18:28,663:iu_debug:	<<<<<< PATCHING IGVFFI3664MJLW To IGVF database with URL https://api.data.igvf.org/IGVFFI3664MJLW and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D0_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1335,
  "md5sum": "69ae1c0c61a899cccd2555aec8e8bdb2",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:28,663:iu_debug:	DRY RUN is enabled.
[4/48] DRY RUN PATCH IGVFFI6014NOEH
2026-10-15 07:18:28,663:iu_debug:	
IN patch()
2026-10-15 07:18:28,663:iu_debug:	>>>>>>GET IGVFFI6014NOEH From DACC with URL https://api.data.igvf.org/IGVFFI6014NOEH/?format=json&frame=edit
2026-10-15 07:18:28,828:iu_debug:	<<<<<< PATCHING IGVFFI6014NOEH To IGVF database with URL https://api.data.igvf.org/IGVFFI6014NOEH and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D0_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1462,
  "md5sum": "0d674cfaf8393f597c3c8c2c8f86096c",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D0_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:28,829:iu_debug:	DRY RUN is enabled.
[5/48] DRY RUN PATCH IGVFFI9935VMFV
2026-10-15 07:18:28,829:iu_debug:	
IN patch()
2026-10-15 07:18:28,829:iu_debug:	>>>>>>GET IGVFFI9935VMFV From DACC with URL https://api.data.igvf.org/IGVFFI9935VMFV/?format=json&frame=edit
2026-10-15 07:18:28,991:iu_debug:	<<<<<< PATCHING IGVFFI9935VMFV To IGVF database with URL https://api.data.igvf.org/IGVFFI9935VMFV and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D1_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1330,
  "md5sum": "baff223b728433f5ccf5f17902fbcdda",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:28,991:iu_debug:	DRY RUN is enabled.
[6/48] DRY RUN PATCH IGVFFI6690AXDH
2026-10-15 07:18:28,991:iu_debug:	
IN patch()
2026-10-15 07:18:28,992:iu_debug:	>>>>>>GET IGVFFI6690AXDH From DACC with URL https://api.data.igvf.org/IGVFFI6690AXDH/?format=json&frame=edit
2026-10-15 07:18:29,287:iu_debug:	<<<<<< PATCHING IGVFFI6690AXDH To IGVF database with URL https://api.data.igvf.org/IGVFFI6690AXDH and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D1_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "853862e522a1e68ced73f135ad8aeb73",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:29,287:iu_debug:	DRY RUN is enabled.
[7/48] DRY RUN PATCH IGVFFI2772BUCG
2026-10-15 07:18:29,287:iu_debug:	
IN patch()
2026-10-15 07:18:29,287:iu_debug:	>>>>>>GET IGVFFI2772BUCG From DACC with URL https://api.data.igvf.org/IGVFFI2772BUCG/?format=json&frame=edit
2026-10-15 07:18:29,442:iu_debug:	<<<<<< PATCHING IGVFFI2772BUCG To IGVF database with URL https://api.data.igvf.org/IGVFFI2772BUCG and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D1_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "0d49bbe9cc466787fcee5cf2b48a2b9f",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:29,442:iu_debug:	DRY RUN is enabled.
[8/48] DRY RUN PATCH IGVFFI0810QEPN
2026-10-15 07:18:29,443:iu_debug:	
IN patch()
2026-10-15 07:18:29,443:iu_debug:	>>>>>>GET IGVFFI0810QEPN From DACC with URL https://api.data.igvf.org/IGVFFI0810QEPN/?format=json&frame=edit
2026-10-15 07:18:29,609:iu_debug:	<<<<<< PATCHING IGVFFI0810QEPN To IGVF database with URL https://api.data.igvf.org/IGVFFI0810QEPN and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D1_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "17b91b1dc303fdebef453401e9f47ce0",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D1_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:29,610:iu_debug:	DRY RUN is enabled.
[9/48] DRY RUN PATCH IGVFFI8048DGMK
2026-10-15 07:18:29,610:iu_debug:	
IN patch()
2026-10-15 07:18:29,610:iu_debug:	>>>>>>GET IGVFFI8048DGMK From DACC with URL https://api.data.igvf.org/IGVFFI8048DGMK/?format=json&frame=edit
2026-10-15 07:18:29,768:iu_debug:	<<<<<< PATCHING IGVFFI8048DGMK To IGVF database with URL https://api.data.igvf.org/IGVFFI8048DGMK and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D2_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1331,
  "md5sum": "d010a6e905ea42f8b9447b18f3c267a1",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:29,768:iu_debug:	DRY RUN is enabled.
[10/48] DRY RUN PATCH IGVFFI8873GOQM
2026-10-15 07:18:29,768:iu_debug:	
IN patch()
2026-10-15 07:18:29,768:iu_debug:	>>>>>>GET IGVFFI8873GOQM From DACC with URL https://api.data.igvf.org/IGVFFI8873GOQM/?format=json&frame=edit
2026-10-15 07:18:29,926:iu_debug:	<<<<<< PATCHING IGVFFI8873GOQM To IGVF database with URL https://api.data.igvf.org/IGVFFI8873GOQM and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D2_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1463,
  "md5sum": "12af6eada8f23cf658cf577d0d1dee42",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:29,927:iu_debug:	DRY RUN is enabled.
[11/48] DRY RUN PATCH IGVFFI4312QTKI
2026-10-15 07:18:29,927:iu_debug:	
IN patch()
2026-10-15 07:18:29,927:iu_debug:	>>>>>>GET IGVFFI4312QTKI From DACC with URL https://api.data.igvf.org/IGVFFI4312QTKI/?format=json&frame=edit
2026-10-15 07:18:30,095:iu_debug:	<<<<<< PATCHING IGVFFI4312QTKI To IGVF database with URL https://api.data.igvf.org/IGVFFI4312QTKI and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D2_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "a4175d31fa459ad98231f40e321749dd",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:30,095:iu_debug:	DRY RUN is enabled.
[12/48] DRY RUN PATCH IGVFFI9788HFUZ
2026-10-15 07:18:30,095:iu_debug:	
IN patch()
2026-10-15 07:18:30,095:iu_debug:	>>>>>>GET IGVFFI9788HFUZ From DACC with URL https://api.data.igvf.org/IGVFFI9788HFUZ/?format=json&frame=edit
2026-10-15 07:18:30,250:iu_debug:	<<<<<< PATCHING IGVFFI9788HFUZ To IGVF database with URL https://api.data.igvf.org/IGVFFI9788HFUZ and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D2_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1463,
  "md5sum": "a01a15439b9e57e728af8671b2891152",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D2_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:30,251:iu_debug:	DRY RUN is enabled.
[13/48] DRY RUN PATCH IGVFFI7854RUAJ
2026-10-15 07:18:30,251:iu_debug:	
IN patch()
2026-10-15 07:18:30,251:iu_debug:	>>>>>>GET IGVFFI7854RUAJ From DACC with URL https://api.data.igvf.org/IGVFFI7854RUAJ/?format=json&frame=edit
2026-10-15 07:18:30,416:iu_debug:	<<<<<< PATCHING IGVFFI7854RUAJ To IGVF database with URL https://api.data.igvf.org/IGVFFI7854RUAJ and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D3_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1334,
  "md5sum": "5cba90572f891f3d162b67508b710c61",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:30,416:iu_debug:	DRY RUN is enabled.
[14/48] DRY RUN PATCH IGVFFI9066UFCN
2026-10-15 07:18:30,416:iu_debug:	
IN patch()
2026-10-15 07:18:30,416:iu_debug:	>>>>>>GET IGVFFI9066UFCN From DACC with URL https://api.data.igvf.org/IGVFFI9066UFCN/?format=json&frame=edit
2026-10-15 07:18:30,578:iu_debug:	<<<<<< PATCHING IGVFFI9066UFCN To IGVF database with URL https://api.data.igvf.org/IGVFFI9066UFCN and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D3_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "8ee334255160173fa52ebdb60cd5f50c",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:30,579:iu_debug:	DRY RUN is enabled.
[15/48] DRY RUN PATCH IGVFFI2187SJIK
2026-10-15 07:18:30,579:iu_debug:	
IN patch()
2026-10-15 07:18:30,579:iu_debug:	>>>>>>GET IGVFFI2187SJIK From DACC with URL https://api.data.igvf.org/IGVFFI2187SJIK/?format=json&frame=edit
2026-10-15 07:18:30,735:iu_debug:	<<<<<< PATCHING IGVFFI2187SJIK To IGVF database with URL https://api.data.igvf.org/IGVFFI2187SJIK and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary1_D3_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "9ed92765c25c4fb3ffe85be192acc5d5",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:30,735:iu_debug:	DRY RUN is enabled.
[16/48] DRY RUN PATCH IGVFFI8158UWTN
2026-10-15 07:18:30,735:iu_debug:	
IN patch()
2026-10-15 07:18:30,735:iu_debug:	>>>>>>GET IGVFFI8158UWTN From DACC with URL https://api.data.igvf.org/IGVFFI8158UWTN/?format=json&frame=edit
2026-10-15 07:18:30,891:iu_debug:	<<<<<< PATCHING IGVFFI8158UWTN To IGVF database with URL https://api.data.igvf.org/IGVFFI8158UWTN and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary1_D3_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1462,
  "md5sum": "c5e12968b521167016914f778f9ced3d",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary1_D3_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:30,892:iu_debug:	DRY RUN is enabled.
[17/48] DRY RUN PATCH IGVFFI4622KMQE
2026-10-15 07:18:30,892:iu_debug:	
IN patch()
2026-10-15 07:18:30,892:iu_debug:	>>>>>>GET IGVFFI4622KMQE From DACC with URL https://api.data.igvf.org/IGVFFI4622KMQE/?format=json&frame=edit
2026-10-15 07:18:31,051:iu_debug:	<<<<<< PATCHING IGVFFI4622KMQE To IGVF database with URL https://api.data.igvf.org/IGVFFI4622KMQE and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D0_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1334,
  "md5sum": "9630cfd0adc73309449af037be888a3a",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:31,051:iu_debug:	DRY RUN is enabled.
[18/48] DRY RUN PATCH IGVFFI8076FSQO
2026-10-15 07:18:31,051:iu_debug:	
IN patch()
2026-10-15 07:18:31,051:iu_debug:	>>>>>>GET IGVFFI8076FSQO From DACC with URL https://api.data.igvf.org/IGVFFI8076FSQO/?format=json&frame=edit
2026-10-15 07:18:31,206:iu_debug:	<<<<<< PATCHING IGVFFI8076FSQO To IGVF database with URL https://api.data.igvf.org/IGVFFI8076FSQO and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D0_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1468,
  "md5sum": "e734696635c74ad4909e17357c4f50e3",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:31,207:iu_debug:	DRY RUN is enabled.
[19/48] DRY RUN PATCH IGVFFI2175EAGA
2026-10-15 07:18:31,207:iu_debug:	
IN patch()
2026-10-15 07:18:31,207:iu_debug:	>>>>>>GET IGVFFI2175EAGA From DACC with URL https://api.data.igvf.org/IGVFFI2175EAGA/?format=json&frame=edit
2026-10-15 07:18:31,361:iu_debug:	<<<<<< PATCHING IGVFFI2175EAGA To IGVF database with URL https://api.data.igvf.org/IGVFFI2175EAGA and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D0_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1334,
  "md5sum": "faa2095eca4dee01093ccdb7ad4d2cdd",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:31,362:iu_debug:	DRY RUN is enabled.
[20/48] DRY RUN PATCH IGVFFI2434RZUD
2026-10-15 07:18:31,362:iu_debug:	
IN patch()
2026-10-15 07:18:31,362:iu_debug:	>>>>>>GET IGVFFI2434RZUD From DACC with URL https://api.data.igvf.org/IGVFFI2434RZUD/?format=json&frame=edit
2026-10-15 07:18:31,512:iu_debug:	<<<<<< PATCHING IGVFFI2434RZUD To IGVF database with URL https://api.data.igvf.org/IGVFFI2434RZUD and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D0_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1466,
  "md5sum": "dbb45d27583479e23b284e21f52997b8",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D0_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:31,512:iu_debug:	DRY RUN is enabled.
[21/48] DRY RUN PATCH IGVFFI8770RPJT
2026-10-15 07:18:31,513:iu_debug:	
IN patch()
2026-10-15 07:18:31,513:iu_debug:	>>>>>>GET IGVFFI8770RPJT From DACC with URL https://api.data.igvf.org/IGVFFI8770RPJT/?format=json&frame=edit
2026-10-15 07:18:31,669:iu_debug:	<<<<<< PATCHING IGVFFI8770RPJT To IGVF database with URL https://api.data.igvf.org/IGVFFI8770RPJT and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D1_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "1960b0b7624bd5421d58cb83fead87aa",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:31,670:iu_debug:	DRY RUN is enabled.
[22/48] DRY RUN PATCH IGVFFI1488HGCB
2026-10-15 07:18:31,670:iu_debug:	
IN patch()
2026-10-15 07:18:31,670:iu_debug:	>>>>>>GET IGVFFI1488HGCB From DACC with URL https://api.data.igvf.org/IGVFFI1488HGCB/?format=json&frame=edit
2026-10-15 07:18:31,826:iu_debug:	<<<<<< PATCHING IGVFFI1488HGCB To IGVF database with URL https://api.data.igvf.org/IGVFFI1488HGCB and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D1_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "440325f9ea678587ca6b96ceea1323cc",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:31,826:iu_debug:	DRY RUN is enabled.
[23/48] DRY RUN PATCH IGVFFI3154UJVP
2026-10-15 07:18:31,826:iu_debug:	
IN patch()
2026-10-15 07:18:31,826:iu_debug:	>>>>>>GET IGVFFI3154UJVP From DACC with URL https://api.data.igvf.org/IGVFFI3154UJVP/?format=json&frame=edit
2026-10-15 07:18:31,992:iu_debug:	<<<<<< PATCHING IGVFFI3154UJVP To IGVF database with URL https://api.data.igvf.org/IGVFFI3154UJVP and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D1_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1331,
  "md5sum": "fb8ce0d2d54727116678f5f20bc6d56a",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:31,992:iu_debug:	DRY RUN is enabled.
[24/48] DRY RUN PATCH IGVFFI2812VJVV
2026-10-15 07:18:31,992:iu_debug:	
IN patch()
2026-10-15 07:18:31,992:iu_debug:	>>>>>>GET IGVFFI2812VJVV From DACC with URL https://api.data.igvf.org/IGVFFI2812VJVV/?format=json&frame=edit
2026-10-15 07:18:32,146:iu_debug:	<<<<<< PATCHING IGVFFI2812VJVV To IGVF database with URL https://api.data.igvf.org/IGVFFI2812VJVV and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D1_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "2a5a0e3fa302cea70cc4e88205720025",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D1_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:32,146:iu_debug:	DRY RUN is enabled.
[25/48] DRY RUN PATCH IGVFFI2347AMLW
2026-10-15 07:18:32,147:iu_debug:	
IN patch()
2026-10-15 07:18:32,147:iu_debug:	>>>>>>GET IGVFFI2347AMLW From DACC with URL https://api.data.igvf.org/IGVFFI2347AMLW/?format=json&frame=edit
2026-10-15 07:18:32,309:iu_debug:	<<<<<< PATCHING IGVFFI2347AMLW To IGVF database with URL https://api.data.igvf.org/IGVFFI2347AMLW and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D2_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1332,
  "md5sum": "48b96d04630d2f254a37d33bd95c5443",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:32,309:iu_debug:	DRY RUN is enabled.
[26/48] DRY RUN PATCH IGVFFI7718VJVV
2026-10-15 07:18:32,309:iu_debug:	
IN patch()
2026-10-15 07:18:32,309:iu_debug:	>>>>>>GET IGVFFI7718VJVV From DACC with URL https://api.data.igvf.org/IGVFFI7718VJVV/?format=json&frame=edit
2026-10-15 07:18:32,464:iu_debug:	<<<<<< PATCHING IGVFFI7718VJVV To IGVF database with URL https://api.data.igvf.org/IGVFFI7718VJVV and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D2_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "5f020bf5d52145a225aaab97320e77ba",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:32,464:iu_debug:	DRY RUN is enabled.
[27/48] DRY RUN PATCH IGVFFI5545EVXX
2026-10-15 07:18:32,464:iu_debug:	
IN patch()
2026-10-15 07:18:32,464:iu_debug:	>>>>>>GET IGVFFI5545EVXX From DACC with URL https://api.data.igvf.org/IGVFFI5545EVXX/?format=json&frame=edit
2026-10-15 07:18:32,619:iu_debug:	<<<<<< PATCHING IGVFFI5545EVXX To IGVF database with URL https://api.data.igvf.org/IGVFFI5545EVXX and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D2_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1330,
  "md5sum": "875527aec4acdbab8b2e766b82c69941",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:32,620:iu_debug:	DRY RUN is enabled.
[28/48] DRY RUN PATCH IGVFFI0654IPLV
2026-10-15 07:18:32,620:iu_debug:	
IN patch()
2026-10-15 07:18:32,620:iu_debug:	>>>>>>GET IGVFFI0654IPLV From DACC with URL https://api.data.igvf.org/IGVFFI0654IPLV/?format=json&frame=edit
2026-10-15 07:18:32,773:iu_debug:	<<<<<< PATCHING IGVFFI0654IPLV To IGVF database with URL https://api.data.igvf.org/IGVFFI0654IPLV and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D2_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "73bb583cc2e7d5d5fb4dcd458731e785",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D2_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:32,774:iu_debug:	DRY RUN is enabled.
[29/48] DRY RUN PATCH IGVFFI8139DPYA
2026-10-15 07:18:32,774:iu_debug:	
IN patch()
2026-10-15 07:18:32,774:iu_debug:	>>>>>>GET IGVFFI8139DPYA From DACC with URL https://api.data.igvf.org/IGVFFI8139DPYA/?format=json&frame=edit
2026-10-15 07:18:32,936:iu_debug:	<<<<<< PATCHING IGVFFI8139DPYA To IGVF database with URL https://api.data.igvf.org/IGVFFI8139DPYA and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D3_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1334,
  "md5sum": "80d1cb6a893357bb5c66976212640c07",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:32,936:iu_debug:	DRY RUN is enabled.
[30/48] DRY RUN PATCH IGVFFI1456NTTN
2026-10-15 07:18:32,937:iu_debug:	
IN patch()
2026-10-15 07:18:32,937:iu_debug:	>>>>>>GET IGVFFI1456NTTN From DACC with URL https://api.data.igvf.org/IGVFFI1456NTTN/?format=json&frame=edit
2026-10-15 07:18:33,097:iu_debug:	<<<<<< PATCHING IGVFFI1456NTTN To IGVF database with URL https://api.data.igvf.org/IGVFFI1456NTTN and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D3_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "40bdd27e60fca03f7b19b06922593b5f",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:33,097:iu_debug:	DRY RUN is enabled.
[31/48] DRY RUN PATCH IGVFFI8337KFXH
2026-10-15 07:18:33,097:iu_debug:	
IN patch()
2026-10-15 07:18:33,097:iu_debug:	>>>>>>GET IGVFFI8337KFXH From DACC with URL https://api.data.igvf.org/IGVFFI8337KFXH/?format=json&frame=edit
2026-10-15 07:18:33,274:iu_debug:	<<<<<< PATCHING IGVFFI8337KFXH To IGVF database with URL https://api.data.igvf.org/IGVFFI8337KFXH and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary2_D3_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1335,
  "md5sum": "0b1f2a172cbae53ff370d158b1e7804b",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:33,275:iu_debug:	DRY RUN is enabled.
[32/48] DRY RUN PATCH IGVFFI6865SJSV
2026-10-15 07:18:33,275:iu_debug:	
IN patch()
2026-10-15 07:18:33,275:iu_debug:	>>>>>>GET IGVFFI6865SJSV From DACC with URL https://api.data.igvf.org/IGVFFI6865SJSV/?format=json&frame=edit
2026-10-15 07:18:33,438:iu_debug:	<<<<<< PATCHING IGVFFI6865SJSV To IGVF database with URL https://api.data.igvf.org/IGVFFI6865SJSV and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary2_D3_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1465,
  "md5sum": "fa9c53d77adacf3150814cee0c18c1da",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary2_D3_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:33,438:iu_debug:	DRY RUN is enabled.
[33/48] DRY RUN PATCH IGVFFI5741OOFJ
2026-10-15 07:18:33,438:iu_debug:	
IN patch()
2026-10-15 07:18:33,438:iu_debug:	>>>>>>GET IGVFFI5741OOFJ From DACC with URL https://api.data.igvf.org/IGVFFI5741OOFJ/?format=json&frame=edit
2026-10-15 07:18:33,593:iu_debug:	<<<<<< PATCHING IGVFFI5741OOFJ To IGVF database with URL https://api.data.igvf.org/IGVFFI5741OOFJ and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D0_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1332,
  "md5sum": "f736e180dfbf4e1b46e8db3e580a725c",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:33,593:iu_debug:	DRY RUN is enabled.
[34/48] DRY RUN PATCH IGVFFI3023CUFB
2026-10-15 07:18:33,593:iu_debug:	
IN patch()
2026-10-15 07:18:33,593:iu_debug:	>>>>>>GET IGVFFI3023CUFB From DACC with URL https://api.data.igvf.org/IGVFFI3023CUFB/?format=json&frame=edit
2026-10-15 07:18:33,753:iu_debug:	<<<<<< PATCHING IGVFFI3023CUFB To IGVF database with URL https://api.data.igvf.org/IGVFFI3023CUFB and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D0_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1461,
  "md5sum": "0b686b82e22dc8746620ab4a75fdeb6c",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:33,754:iu_debug:	DRY RUN is enabled.
[35/48] DRY RUN PATCH IGVFFI9717GCXG
2026-10-15 07:18:33,754:iu_debug:	
IN patch()
2026-10-15 07:18:33,754:iu_debug:	>>>>>>GET IGVFFI9717GCXG From DACC with URL https://api.data.igvf.org/IGVFFI9717GCXG/?format=json&frame=edit
2026-10-15 07:18:33,904:iu_debug:	<<<<<< PATCHING IGVFFI9717GCXG To IGVF database with URL https://api.data.igvf.org/IGVFFI9717GCXG and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D0_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "d44cc9060b28a8821b68f868fec81775",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:33,904:iu_debug:	DRY RUN is enabled.
[36/48] DRY RUN PATCH IGVFFI5837FVFN
2026-10-15 07:18:33,904:iu_debug:	
IN patch()
2026-10-15 07:18:33,904:iu_debug:	>>>>>>GET IGVFFI5837FVFN From DACC with URL https://api.data.igvf.org/IGVFFI5837FVFN/?format=json&frame=edit
2026-10-15 07:18:34,062:iu_debug:	<<<<<< PATCHING IGVFFI5837FVFN To IGVF database with URL https://api.data.igvf.org/IGVFFI5837FVFN and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D0_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1463,
  "md5sum": "d03ff95ed56b7421ceca192ea67f05e3",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D0_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:34,062:iu_debug:	DRY RUN is enabled.
[37/48] DRY RUN PATCH IGVFFI6943RGAA
2026-10-15 07:18:34,062:iu_debug:	
IN patch()
2026-10-15 07:18:34,063:iu_debug:	>>>>>>GET IGVFFI6943RGAA From DACC with URL https://api.data.igvf.org/IGVFFI6943RGAA/?format=json&frame=edit
2026-10-15 07:18:34,220:iu_debug:	<<<<<< PATCHING IGVFFI6943RGAA To IGVF database with URL https://api.data.igvf.org/IGVFFI6943RGAA and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D1_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1336,
  "md5sum": "09b56465120c2a1326a788c68dcd2611",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:34,221:iu_debug:	DRY RUN is enabled.
[38/48] DRY RUN PATCH IGVFFI6783SCTB
2026-10-15 07:18:34,221:iu_debug:	
IN patch()
2026-10-15 07:18:34,221:iu_debug:	>>>>>>GET IGVFFI6783SCTB From DACC with URL https://api.data.igvf.org/IGVFFI6783SCTB/?format=json&frame=edit
2026-10-15 07:18:34,374:iu_debug:	<<<<<< PATCHING IGVFFI6783SCTB To IGVF database with URL https://api.data.igvf.org/IGVFFI6783SCTB and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D1_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1465,
  "md5sum": "4f43cb50bb8d712e0497941ff60b1665",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:34,374:iu_debug:	DRY RUN is enabled.
[39/48] DRY RUN PATCH IGVFFI3958VUUO
2026-10-15 07:18:34,374:iu_debug:	
IN patch()
2026-10-15 07:18:34,374:iu_debug:	>>>>>>GET IGVFFI3958VUUO From DACC with URL https://api.data.igvf.org/IGVFFI3958VUUO/?format=json&frame=edit
2026-10-15 07:18:34,528:iu_debug:	<<<<<< PATCHING IGVFFI3958VUUO To IGVF database with URL https://api.data.igvf.org/IGVFFI3958VUUO and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D1_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1334,
  "md5sum": "5af5c6c2587e25772dd6d16190f771a0",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:34,529:iu_debug:	DRY RUN is enabled.
[40/48] DRY RUN PATCH IGVFFI4268ZJMH
2026-10-15 07:18:34,529:iu_debug:	
IN patch()
2026-10-15 07:18:34,529:iu_debug:	>>>>>>GET IGVFFI4268ZJMH From DACC with URL https://api.data.igvf.org/IGVFFI4268ZJMH/?format=json&frame=edit
2026-10-15 07:18:34,677:iu_debug:	<<<<<< PATCHING IGVFFI4268ZJMH To IGVF database with URL https://api.data.igvf.org/IGVFFI4268ZJMH and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D1_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1465,
  "md5sum": "0ce6450f4b79bc8d2ff81f28cbb97ff2",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D1_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:34,677:iu_debug:	DRY RUN is enabled.
[41/48] DRY RUN PATCH IGVFFI4966LKOU
2026-10-15 07:18:34,677:iu_debug:	
IN patch()
2026-10-15 07:18:34,677:iu_debug:	>>>>>>GET IGVFFI4966LKOU From DACC with URL https://api.data.igvf.org/IGVFFI4966LKOU/?format=json&frame=edit
2026-10-15 07:18:34,836:iu_debug:	<<<<<< PATCHING IGVFFI4966LKOU To IGVF database with URL https://api.data.igvf.org/IGVFFI4966LKOU and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D2_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "539978339e061580d592f5eccdf25b78",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:34,836:iu_debug:	DRY RUN is enabled.
[42/48] DRY RUN PATCH IGVFFI9383WWAH
2026-10-15 07:18:34,836:iu_debug:	
IN patch()
2026-10-15 07:18:34,836:iu_debug:	>>>>>>GET IGVFFI9383WWAH From DACC with URL https://api.data.igvf.org/IGVFFI9383WWAH/?format=json&frame=edit
2026-10-15 07:18:34,991:iu_debug:	<<<<<< PATCHING IGVFFI9383WWAH To IGVF database with URL https://api.data.igvf.org/IGVFFI9383WWAH and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D2_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1465,
  "md5sum": "c85e91e34612113fd2f5b3aca774a435",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:34,991:iu_debug:	DRY RUN is enabled.
[43/48] DRY RUN PATCH IGVFFI5241GSQF
2026-10-15 07:18:34,991:iu_debug:	
IN patch()
2026-10-15 07:18:34,991:iu_debug:	>>>>>>GET IGVFFI5241GSQF From DACC with URL https://api.data.igvf.org/IGVFFI5241GSQF/?format=json&frame=edit
2026-10-15 07:18:35,147:iu_debug:	<<<<<< PATCHING IGVFFI5241GSQF To IGVF database with URL https://api.data.igvf.org/IGVFFI5241GSQF and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D2_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1336,
  "md5sum": "d92c87068ca38697bae7bc1b63b20177",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:35,147:iu_debug:	DRY RUN is enabled.
[44/48] DRY RUN PATCH IGVFFI8680AUTL
2026-10-15 07:18:35,147:iu_debug:	
IN patch()
2026-10-15 07:18:35,147:iu_debug:	>>>>>>GET IGVFFI8680AUTL From DACC with URL https://api.data.igvf.org/IGVFFI8680AUTL/?format=json&frame=edit
2026-10-15 07:18:35,304:iu_debug:	<<<<<< PATCHING IGVFFI8680AUTL To IGVF database with URL https://api.data.igvf.org/IGVFFI8680AUTL and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D2_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "a685ab81f917d521995d3f427136dd8a",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D2_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:35,305:iu_debug:	DRY RUN is enabled.
[45/48] DRY RUN PATCH IGVFFI3799GGVF
2026-10-15 07:18:35,305:iu_debug:	
IN patch()
2026-10-15 07:18:35,305:iu_debug:	>>>>>>GET IGVFFI3799GGVF From DACC with URL https://api.data.igvf.org/IGVFFI3799GGVF/?format=json&frame=edit
2026-10-15 07:18:35,465:iu_debug:	<<<<<< PATCHING IGVFFI3799GGVF To IGVF database with URL https://api.data.igvf.org/IGVFFI3799GGVF and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D3_biorep1; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "dfec8fb9ca8162e4f97da8601ea8c31b",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep1_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:35,465:iu_debug:	DRY RUN is enabled.
[46/48] DRY RUN PATCH IGVFFI8059GUMK
2026-10-15 07:18:35,466:iu_debug:	
IN patch()
2026-10-15 07:18:35,466:iu_debug:	>>>>>>GET IGVFFI8059GUMK From DACC with URL https://api.data.igvf.org/IGVFFI8059GUMK/?format=json&frame=edit
2026-10-15 07:18:35,650:iu_debug:	<<<<<< PATCHING IGVFFI8059GUMK To IGVF database with URL https://api.data.igvf.org/IGVFFI8059GUMK and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D3_biorep1; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "6fd07e5b66f8498413827175e8f962d9",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep1_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:35,650:iu_debug:	DRY RUN is enabled.
[47/48] DRY RUN PATCH IGVFFI5429IGZN
2026-10-15 07:18:35,651:iu_debug:	
IN patch()
2026-10-15 07:18:35,651:iu_debug:	>>>>>>GET IGVFFI5429IGZN From DACC with URL https://api.data.igvf.org/IGVFFI5429IGZN/?format=json&frame=edit
2026-10-15 07:18:35,813:iu_debug:	<<<<<< PATCHING IGVFFI5429IGZN To IGVF database with URL https://api.data.igvf.org/IGVFFI5429IGZN and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq GEX seqspec for sublibrary3_D3_biorep2; regenerated using the corrected GEX structural template represented by IGVFFI6906STSO, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1333,
  "md5sum": "426287b333879d78c92be3136d2eff70",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep2_GEX_seqspec.yaml.gz"
}


2026-10-15 07:18:35,813:iu_debug:	DRY RUN is enabled.
[48/48] DRY RUN PATCH IGVFFI2335DNRR
2026-10-15 07:18:35,813:iu_debug:	
IN patch()
2026-10-15 07:18:35,813:iu_debug:	>>>>>>GET IGVFFI2335DNRR From DACC with URL https://api.data.igvf.org/IGVFFI2335DNRR/?format=json&frame=edit
2026-10-15 07:18:35,964:iu_debug:	<<<<<< PATCHING IGVFFI2335DNRR To IGVF database with URL https://api.data.igvf.org/IGVFFI2335DNRR and this payload:

{
  "description": "WTC11-EC 300-gene CC-Perturb-seq gRNA seqspec for sublibrary3_D3_biorep2; regenerated using the corrected gRNA structural template represented by IGVFFI4006RNQZ, while preserving the original sample-specific read/file metadata and date.",
  "file_size": 1464,
  "md5sum": "3c8bd80fb492e4db81e6862d40bf303d",
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz/WTC11-EC_300gene_CC-Perturb-seq_sublibrary3_D3_biorep2_gRNA_seqspec.yaml.gz"
}


2026-10-15 07:18:35,964:iu_debug:	DRY RUN is enabled.

All metadata dry runs passed.
Dry-run marker: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/_work/METADATA_DRY_RUN_PASSED.json

```

## Markdown cell 22

## 11. Replace the existing 48 files in production

## Code cell 23

```python
PRODUCTION_CONFIRMATION = (
    "REPLACE_48_300GENE_SEQSPEC_FILES_IN_PLACE"
)
```

## Code cell 24

```python
def find_aws_executable() -> Path:
    candidates = [
        Path(shutil.which("aws")) if shutil.which("aws") else None,
        IGVF_PYTHON.parent / "aws",
        Path.home() / ".local/bin/aws",
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "AWS CLI was not found in PATH or the igvf_utils environment"
    )


def issue_upload_credentials(
    session: requests.Session,
    accession: str,
) -> dict:
    response = session.post(
        f"{API_BASE}/files/{accession}/@@upload",
        json={},
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        timeout=REQUEST_TIMEOUT_SECONDS,
    )
    if not response.ok:
        raise RuntimeError(
            f"Could not issue upload credentials for {accession}: "
            f"HTTP {response.status_code}\n{response.text[:4000]}"
        )
    payload = response.json()
    credentials = payload["@graph"][0]["upload_credentials"]
    required = {"access_key", "secret_key", "session_token", "upload_url"}
    if not required <= set(credentials):
        raise RuntimeError(
            f"Upload credentials for {accession} are missing keys: "
            f"{sorted(required - set(credentials))}"
        )
    return credentials


def aws_environment(credentials: dict) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update(
        {
            "AWS_ACCESS_KEY_ID": credentials["access_key"],
            "AWS_SECRET_ACCESS_KEY": credentials["secret_key"],
            "AWS_SESSION_TOKEN": credentials["session_token"],
        }
    )
    return environment


def upload_to_s3(
    aws_executable: Path,
    local_path: Path,
    credentials: dict,
    log_path: Path,
) -> dict:
    upload_url = credentials["upload_url"]
    environment = aws_environment(credentials)

    upload_result = subprocess.run(
        [str(aws_executable), "s3", "cp", str(local_path), upload_url, "--only-show-errors"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
        check=False,
    )
    log_path.write_text(upload_result.stdout, encoding="utf-8")
    if upload_result.returncode != 0:
        raise RuntimeError(
            f"S3 upload failed for {local_path}:\n{upload_result.stdout}"
        )

    parsed = urlparse(upload_url)
    bucket = parsed.netloc
    key = parsed.path.lstrip("/")
    head_result = subprocess.run(
        [
            str(aws_executable),
            "s3api",
            "head-object",
            "--bucket",
            bucket,
            "--key",
            key,
            "--output",
            "json",
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
        check=False,
    )

    if head_result.returncode == 0:
        head = json.loads(head_result.stdout)
        head_success = True
        head_error = ""
        content_length = int(head.get("ContentLength", -1))
        etag = str(head.get("ETag", "")).strip('"')
    else:
        # Some upload credentials allow PutObject but not HeadObject. The AWS
        # upload command itself is still authoritative for transfer success.
        head_success = False
        head_error = head_result.stdout.strip()
        content_length = -1
        etag = ""

    return {
        "upload_url": upload_url,
        "s3_head_success": head_success,
        "s3_head_error": head_error,
        "s3_content_length": content_length,
        "s3_etag": etag,
        "upload_log": str(log_path),
    }


def portal_download_md5(
    session: requests.Session,
    record: dict,
    accession: str,
) -> tuple[str, str]:
    at_id = str(record.get("@id", ""))
    if not at_id:
        return "", "record has no @id"
    url = (
        API_BASE.rstrip("/")
        + at_id
        + "@@download/"
        + accession
        + ".yaml.gz"
    )
    try:
        response = session.get(
            url,
            timeout=REQUEST_TIMEOUT_SECONDS,
            allow_redirects=True,
        )
        if not response.ok:
            return "", f"HTTP {response.status_code}"
        return md5_bytes(response.content), ""
    except Exception as error:
        return "", f"{type(error).__name__}: {error}"


if PRODUCTION_CONFIRMATION != REQUIRED_PRODUCTION_CONFIRMATION:
    print("Production replacement is blocked.")
    print("Set PRODUCTION_CONFIRMATION to:")
    print(REQUIRED_PRODUCTION_CONFIRMATION)
else:
    if not DRY_RUN_MARKER.is_file():
        raise RuntimeError("Metadata dry-run marker is missing")

    marker = json.loads(DRY_RUN_MARKER.read_text(encoding="utf-8"))
    if marker["generation_manifest_sha256"] != sha256_file(generation_manifest_path):
        raise RuntimeError("Generation manifest changed after dry run")
    if marker["portal_plan_sha256"] != sha256_file(portal_plan_path):
        raise RuntimeError("Portal update plan changed after dry run")
    if marker["operations_sha256"] != sha256_file(operations_path):
        raise RuntimeError("Update operations changed after dry run")

    for row in generation_manifest.to_dict("records"):
        if marker["files"][row["accession"]] != row["corrected_sha256"]:
            raise RuntimeError(
                f"Corrected file changed after dry run: {row['accession']}"
            )

    aws_executable = find_aws_executable()
    production_results_path = REPORT_DIR / "production_upload_results.tsv"

    previous_rows = []
    if production_results_path.is_file():
        previous_rows = pd.read_csv(
            production_results_path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        ).to_dict("records")
    previous_by_accession = {
        row["accession"]: row
        for row in previous_rows
        if row.get("verification") == "PASS"
    }

    results: list[dict] = []

    for index, operation in enumerate(update_operations, start=1):
        accession = operation["accession"]
        corrected_path = Path(operation["corrected_gz_path"])
        expected_md5 = operation["corrected_md5"]
        expected_size = int(operation["corrected_size"])

        current_record = get_portal_record(
            portal_session,
            accession,
            frame="object",
            database=True,
        )
        current_md5 = str(current_record.get("md5sum", ""))
        current_size = int(current_record.get("file_size", -1) or -1)

        if (
            accession in previous_by_accession
            and current_md5 == expected_md5
            and current_size == expected_size
        ):
            results.append(previous_by_accession[accession])
            print(f"[{index}/48] SKIP verified {accession}", flush=True)
            continue

        if current_md5 not in {
            expected_md5,
            map_by_accession[accession]["md5sum"],
        }:
            raise RuntimeError(
                f"Unexpected current Portal MD5 for {accession}: {current_md5}"
            )

        print(f"[{index}/48] replacing {accession}", flush=True)

        # Request credentials before mutating file metadata. A credential error
        # therefore cannot leave a partially patched record.
        credentials = issue_upload_credentials(portal_session, accession)

        patch_payload = dict(operation["payload"])
        patch_payload.pop("record_id", None)
        patch_response = patch_portal_record(
            portal_session,
            accession,
            patch_payload,
        )

        upload_log = LOG_DIR / f"{accession}_aws_upload.log"
        s3_result = upload_to_s3(
            aws_executable,
            corrected_path,
            credentials,
            upload_log,
        )

        if (
            s3_result["s3_head_success"]
            and s3_result["s3_content_length"] != expected_size
        ):
            raise RuntimeError(
                f"S3 size mismatch for {accession}: "
                f"expected {expected_size}, observed {s3_result['s3_content_length']}"
            )

        # ETag is recorded as supporting evidence but is not treated as a
        # universal content hash because S3 does not guarantee MD5 semantics
        # for every encryption/upload configuration.
        etag = s3_result["s3_etag"]
        etag_md5_match = (
            bool(re.fullmatch(r"[0-9a-fA-F]{32}", etag))
            and etag.lower() == expected_md5.lower()
        )

        verified_record = get_portal_record(
            portal_session,
            accession,
            frame="object",
            database=True,
        )
        observed_md5 = str(verified_record.get("md5sum", ""))
        observed_size = int(verified_record.get("file_size", -1) or -1)
        observed_path = str(verified_record.get("submitted_file_name", ""))

        if observed_md5 != expected_md5:
            raise RuntimeError(
                f"Portal MD5 verification failed for {accession}: {observed_md5}"
            )
        if observed_size != expected_size:
            raise RuntimeError(
                f"Portal size verification failed for {accession}: {observed_size}"
            )
        if observed_path != str(corrected_path):
            raise RuntimeError(
                f"Portal submitted_file_name verification failed for {accession}: "
                f"{observed_path}"
            )

        remote_md5 = ""
        remote_error = ""
        if TRY_PORTAL_DOWNLOAD_VERIFICATION:
            remote_md5, remote_error = portal_download_md5(
                portal_session,
                verified_record,
                accession,
            )
            if remote_md5 and remote_md5 != expected_md5:
                raise RuntimeError(
                    f"Portal download MD5 mismatch for {accession}: {remote_md5}"
                )

        response_path = PORTAL_BACKUP_DIR / f"{accession}_patch_response.json"
        response_path.write_text(
            json.dumps(patch_response, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        result_row = {
            "accession": accession,
            "modality": operation["modality"],
            "action": "REPLACED_IN_PLACE",
            "expected_md5": expected_md5,
            "observed_portal_md5": observed_md5,
            "expected_size": expected_size,
            "observed_portal_size": observed_size,
            "s3_head_success": s3_result["s3_head_success"],
            "s3_head_error": s3_result["s3_head_error"],
            "s3_content_length": s3_result["s3_content_length"],
            "s3_etag": etag,
            "s3_etag_md5_match": etag_md5_match,
            "portal_download_md5": remote_md5,
            "portal_download_error": remote_error,
            "upload_status": str(verified_record.get("upload_status", "")),
            "status": str(verified_record.get("status", "")),
            "response_json": str(response_path),
            "upload_log": str(upload_log),
            "verification": "PASS",
        }
        results.append(result_row)

        pd.DataFrame(results).to_csv(
            production_results_path,
            sep="\t",
            index=False,
            lineterminator="\n",
        )

    production_results = pd.DataFrame(results)
    if len(production_results) != 48 or not production_results["verification"].eq("PASS").all():
        raise RuntimeError(
            f"Production replacement did not finish 48/48. Review {production_results_path}"
        )

    completion_marker = OUTPUT_ROOT / "CORRECTED_SEQSPECS_REPLACED_48_OF_48.txt"
    completion_marker.write_text(
        "All 48 existing ConfigurationFiles were replaced in place with the "
        "full corrected structural templates.\n"
        f"Results: {production_results_path}\n",
        encoding="utf-8",
    )

    print("Production replacement completed 48/48.")
    print("Results:", production_results_path)
    print("Completion marker:", completion_marker)
```

### Text output (stream)

```text
[1/48] replacing IGVFFI9366VIOI
[2/48] replacing IGVFFI5296LCGA
[3/48] replacing IGVFFI3664MJLW
[4/48] replacing IGVFFI6014NOEH
[5/48] replacing IGVFFI9935VMFV
[6/48] replacing IGVFFI6690AXDH
[7/48] replacing IGVFFI2772BUCG
[8/48] replacing IGVFFI0810QEPN
[9/48] replacing IGVFFI8048DGMK
[10/48] replacing IGVFFI8873GOQM
[11/48] replacing IGVFFI4312QTKI
[12/48] replacing IGVFFI9788HFUZ
[13/48] replacing IGVFFI7854RUAJ
[14/48] replacing IGVFFI9066UFCN
[15/48] replacing IGVFFI2187SJIK
[16/48] replacing IGVFFI8158UWTN
[17/48] replacing IGVFFI4622KMQE
[18/48] replacing IGVFFI8076FSQO
[19/48] replacing IGVFFI2175EAGA
[20/48] replacing IGVFFI2434RZUD
[21/48] replacing IGVFFI8770RPJT
[22/48] replacing IGVFFI1488HGCB
[23/48] replacing IGVFFI3154UJVP
[24/48] replacing IGVFFI2812VJVV
[25/48] replacing IGVFFI2347AMLW
[26/48] replacing IGVFFI7718VJVV
[27/48] replacing IGVFFI5545EVXX
[28/48] replacing IGVFFI0654IPLV
[29/48] replacing IGVFFI8139DPYA
[30/48] replacing IGVFFI1456NTTN
[31/48] replacing IGVFFI8337KFXH
[32/48] replacing IGVFFI6865SJSV
[33/48] replacing IGVFFI5741OOFJ
[34/48] replacing IGVFFI3023CUFB
[35/48] replacing IGVFFI9717GCXG
[36/48] replacing IGVFFI5837FVFN
[37/48] replacing IGVFFI6943RGAA
[38/48] replacing IGVFFI6783SCTB
[39/48] replacing IGVFFI3958VUUO
[40/48] replacing IGVFFI4268ZJMH
[41/48] replacing IGVFFI4966LKOU
[42/48] replacing IGVFFI9383WWAH
[43/48] replacing IGVFFI5241GSQF
[44/48] replacing IGVFFI8680AUTL
[45/48] replacing IGVFFI3799GGVF
[46/48] replacing IGVFFI8059GUMK
[47/48] replacing IGVFFI5429IGZN
[48/48] replacing IGVFFI2335DNRR
Production replacement completed 48/48.
Results: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/production_upload_results.tsv
Completion marker: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/CORRECTED_SEQSPECS_REPLACED_48_OF_48.txt

```

## Markdown cell 25

## 12. Final read-only reconciliation

## Code cell 26

```python
final_rows: list[dict] = []

for index, row in enumerate(generation_manifest.to_dict("records"), start=1):
    accession = row["accession"]
    mapping = map_by_accession[accession]
    record = get_portal_record(
        portal_session,
        accession,
        frame="object",
        database=True,
    )

    aliases = record.get("aliases", [])
    if isinstance(aliases, str):
        aliases = [aliases]
    aliases = set(map(str, aliases))

    expected_seqspec_of = {
        value.strip()
        for value in mapping["seqspec_of"].split(",")
        if value.strip()
    }
    observed_seqspec_of = link_accessions(record.get("seqspec_of", []), "IGVFFI")

    expected_file_set_match = re.search(r"IGVFDS[A-Z0-9]+", mapping["file_set"])
    expected_file_set = {expected_file_set_match.group(0)} if expected_file_set_match else set()
    observed_file_set = link_accessions(record.get("file_set", ""), "IGVFDS")

    checks = {
        "md5": str(record.get("md5sum", "")) == row["corrected_md5"],
        "file_size": int(record.get("file_size", -1) or -1) == int(row["corrected_size"]),
        "submitted_file_name": str(record.get("submitted_file_name", "")) == row["corrected_gz_path"],
        "alias": mapping["alias"] in aliases,
        "file_set": observed_file_set == expected_file_set,
        "seqspec_of": observed_seqspec_of == expected_seqspec_of,
        "content_type": str(record.get("content_type", "")) == "seqspec",
        "file_format": str(record.get("file_format", "")) == "yaml",
    }
    failed = sorted(key for key, value in checks.items() if not value)

    final_rows.append(
        {
            "accession": accession,
            "modality": row["modality"],
            "status": str(record.get("status", "")),
            "upload_status": str(record.get("upload_status", "")),
            "corrected_md5": row["corrected_md5"],
            "portal_md5": str(record.get("md5sum", "")),
            "failed_checks": ",".join(failed),
            "verification": "PASS" if not failed else "FAIL",
        }
    )

    if index == 1 or index % 12 == 0 or index == len(generation_manifest):
        print(f"Final reconciliation {index}/{len(generation_manifest)}", flush=True)

final_reconciliation = pd.DataFrame(final_rows)
final_reconciliation_path = REPORT_DIR / "final_portal_reconciliation.tsv"
write_tsv(final_reconciliation, final_reconciliation_path)

print(final_reconciliation.groupby(["modality", "verification"]).size())
print(final_reconciliation["upload_status"].value_counts(dropna=False))
print("Final reconciliation:", final_reconciliation_path)

if final_reconciliation["verification"].eq("PASS").all():
    marker = OUTPUT_ROOT / "PORTAL_METADATA_RECONCILED_48_OF_48.txt"
    marker.write_text(
        "All 48 existing ConfigurationFile records have the corrected MD5, size, "
        "path, alias, file set, and seqspec_of metadata.\n"
        f"Results: {final_reconciliation_path}\n",
        encoding="utf-8",
    )
    print("Metadata completion marker:", marker)
else:
    raise RuntimeError(
        "Final Portal reconciliation failed. Review "
        f"{final_reconciliation_path}"
    )
```

### Text output (stream)

```text
Final reconciliation 1/48
Final reconciliation 12/48
Final reconciliation 24/48
Final reconciliation 36/48
Final reconciliation 48/48
modality  verification
GEX       PASS            24
gRNA      PASS            24
dtype: int64
upload_status
pending    48
Name: count, dtype: int64
Final reconciliation: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/final_portal_reconciliation.tsv
Metadata completion marker: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/PORTAL_METADATA_RECONCILED_48_OF_48.txt

```

## Markdown cell 27

## 13. Output summary

## Code cell 28

```python
outputs = [
    ("original source inventory", REPORT_DIR / "original_source_inventory.tsv"),
    ("generation manifest", REPORT_DIR / "corrected_generation_manifest.tsv"),
    ("structural change report", REPORT_DIR / "approved_structural_change_report.tsv"),
    ("structural validation", REPORT_DIR / "structural_validation.tsv"),
    ("seqspec CLI validation", REPORT_DIR / "seqspec_cli_validation.tsv"),
    ("package manifest", REPORT_DIR / "package_manifest.json"),
    ("portal update plan", REPORT_DIR / "portal_update_plan.tsv"),
    ("production upload results", REPORT_DIR / "production_upload_results.tsv"),
    ("final reconciliation", REPORT_DIR / "final_portal_reconciliation.tsv"),
    ("resubmission package", OUTPUT_ROOT / "corrected_300gene_seqspec_resubmission_package.zip"),
]

summary = pd.DataFrame(
    [
        {
            "output": label,
            "path": str(path),
            "exists": path.is_file(),
            "size_bytes": path.stat().st_size if path.is_file() else 0,
        }
        for label, path in outputs
    ]
)
print(summary.to_string(index=False))
print("\nCorrected YAML files:", YAML_DIR)
print("Corrected YAML.GZ files:", YAML_GZ_DIR)
print("Logs:", LOG_DIR)
```

### Text output (stream)

```text
                   output                                                                                                                                                          path  exists  size_bytes
original source inventory              /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/original_source_inventory.tsv    True       17349
      generation manifest          /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/corrected_generation_manifest.tsv    True       50830
 structural change report      /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/approved_structural_change_report.tsv    True      109904
    structural validation                  /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/structural_validation.tsv    True        1415
   seqspec CLI validation                 /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/seqspec_cli_validation.tsv    True       19724
         package manifest                      /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/package_manifest.json    True         363
       portal update plan                     /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/portal_update_plan.tsv    True       30663
production upload results              /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/production_upload_results.tsv    True       25973
     final reconciliation            /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/reports/final_portal_reconciliation.tsv    True        5444
     resubmission package /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/corrected_300gene_seqspec_resubmission_package.zip    True      195053

Corrected YAML files: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml
Corrected YAML.GZ files: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/yaml_gz
Logs: /oak/stanford/groups/engreitz/Users/tri/300genes_IGVF_demux/seqspec_two_template_exact/corrected_templates/logs

```

## Code cell 29

```python

```


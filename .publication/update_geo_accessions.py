"""Set study GEO accessions without changing calculations or saved outputs."""
import ast
import copy
import hashlib
import json
from pathlib import Path

NOTEBOOKS = {
    'notebooks/main/Figure_03.ipynb': ('GSE352939', 'Small molecule Hippo scRNA-Seq'),
    'notebooks/main/Figure_04.ipynb': ('GSE352939', 'Small molecule Hippo scRNA-Seq'),
    'notebooks/supplementary/Figure_S02.ipynb': ('GSE352893', 'Parse v3 CC-Perturb-seq, iPSC-K562 mix'),
}

def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'Expected one occurrence: {old!r}; found {text.count(old)}')
    return text.replace(old, new, 1)

def snapshot(nb):
    obj = copy.deepcopy(nb)
    for cell in obj['cells']:
        cell['source'] = None
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False).encode()).hexdigest()

reports = []
for name, (accession, title) in NOTEBOOKS.items():
    path = Path(name)
    nb = json.loads(path.read_text())
    original = copy.deepcopy(nb)
    old_hash = snapshot(nb)
    setting_count = 0
    for cell in nb['cells']:
        text = ''.join(cell['source'])
        if cell['cell_type'] == 'code' and 'GEO_ACCESSION = ""' in text:
            text = replace_once(text,
                'GEO_ACCESSION = ""           # A public GSE... or GSM... accession for this experiment',
                f'GEO_ACCESSION = "{accession}"  # {title}')
            cell['source'] = text.splitlines(keepends=True)
            setting_count += 1
        elif cell['cell_type'] == 'markdown' and 'enter its public `GEO_ACCESSION`' in text:
            text = replace_once(text,
                'set `DATA_SOURCE="geo"` and enter its public `GEO_ACCESSION` (`GSE...` or `GSM...`).',
                f'set `DATA_SOURCE="geo"`; `GEO_ACCESSION="{accession}"` is already configured for [{title}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}).')
            cell['source'] = text.splitlines(keepends=True)
    assert setting_count == 1, name
    for before, after in zip(original['cells'], nb['cells']):
        if before['cell_type'] != 'code':
            continue
        left, right = ast.parse(''.join(before['source'])), ast.parse(''.join(after['source']))
        for node in ast.walk(right):
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'GEO_ACCESSION' for t in node.targets):
                assert ast.literal_eval(node.value) == accession
                node.value = ast.Constant(value='')
        assert ast.dump(left, include_attributes=False) == ast.dump(right, include_attributes=False), name
        compile(''.join(after['source']), name, 'exec')
    import nbformat
    nbformat.validate(nbformat.from_dict(nb))
    path.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + '\n')
    assert snapshot(json.loads(path.read_text())) == old_hash, name
    reports.append({'path': name, 'accession': accession, 'default_mode': 'precomputed',
                    'outputs_execution_counts_metadata_sha256': old_hash,
                    'outputs_preserved': sum(len(c.get('outputs', [])) for c in nb['cells'])})

path = Path('README.md')
s = path.read_text()
s = replace_once(s,
    'checks these files. The Hippo and mixed-cell integrated GEO H5MUs belong in\n`data/geo/`; their download URLs are not configured yet.',
    'checks these files. Hippo and mixed-cell panels use compact precomputed tables\nby default. Optional regeneration downloads the submitted integrated H5MUs from\nGEO: **GSE352939** (Hippo) and **GSE352893** (iPSC-K562 mix).')
s = replace_once(s,
    '**Publication readiness:** adding the download code does not populate the\ncompanion bundle or the GEO files. Until the actual companion data is committed\nand the GEO files are available, affected notebooks stop with explicit missing-\ninput errors rather than using an inaccessible institutional fallback.',
    '**Publication readiness:** default plotting requires the companion data and\nprecomputed tables, not the large GEO objects. Optional GEO regeneration requires\npublicly accessible processed supplementary files. An assigned accession alone\ndoes not guarantee that downloads are released. Missing inputs produce explicit\nerrors; no workflow falls back to institutional filesystem paths.')
s = replace_once(s,
    'companion-file checksums and GEO-file availability; it reports incomplete data',
    'companion and precomputed plotting-table checksums; it reports incomplete data')
s = replace_once(s,
    'needed for this route. Set `DATA_SOURCE="geo"` and `GEO_ACCESSION` in the notebook\nto regenerate those tables from the submitted integrated H5MU; use `local_file`\nmode before GEO release. Regenerated data are kept separately in `data/derived/`.',
    '''needed for this route. The optional regeneration settings are already configured:

| Experiment | GEO accession | Notebook / panels |
|---|---|---|
| Parse v3 CC-Perturb-seq, iPSC-K562 mix | [GSE352893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352893) | `Figure_S02.ipynb`, panels i, j, l, m, n |
| Small molecule Hippo scRNA-Seq | [GSE352939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352939) | `Figure_03.ipynb` panel f and `Figure_04.ipynb` panel c |

To download the submitted object and regenerate plotting inputs, change only:

```python
DATA_SOURCE = "geo"
```

The notebook retains the appropriate `GEO_ACCESSION` above. `GEO_FILENAME=None`
automatically discovers the processed supplementary file; set its exact filename
only if multiple matching files are reported. Released records need no GEO login.
An unreleased record requires `DATA_SOURCE="local_file"` and `GEO_LOCAL_FILE`
pointing to the submitted integrated H5MU until the files become public.

Downloads are cached in `data/geo/downloads/`. Regenerated tables go to
`data/derived/geo_plots/` and are compared with published tables when available;
the published precomputed bundle is never overwritten. Regeneration starts from
the submitted processed RNA plus guide/program modalities, not FASTQ files; it
does not refit cNMF, recluster cells, or reprocess sequencing reads.

Command-line equivalents:

```bash
python scripts/prepare_geo_plot_data.py --dataset ipsc_k562_mix --geo-accession GSE352893
python scripts/prepare_geo_plot_data.py --dataset hippo_timecourse --geo-accession GSE352939
```''')
path.write_text(s)

path = Path('docs/PUBLIC_DATA.md')
s = path.read_text()
s = replace_once(s,
    'GEO_ACCESSION = "GSE123456"  # Replace with this experiment\'s public GSE or GSM',
    'GEO_ACCESSION = "GSE352939"  # Hippo: Figures 3f/4c; Figure S2 uses GSE352893')
s = replace_once(s,
    'The example accession above is illustrative, not an accession assigned to this\nstudy. Public study accessions have not yet been supplied in this repository.',
    '''The notebooks already contain the study-specific accessions below; only
`DATA_SOURCE` needs to change from `"precomputed"` to `"geo"`.

| Dataset | Accession | Notebook |
|---|---|---|
| Parse v3 CC-Perturb-seq, iPSC-K562 mix | [GSE352893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352893) | Figure S2 |
| Small molecule Hippo scRNA-Seq | [GSE352939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352939) | Figures 3 and 4 |

Accessions identify the submissions; anonymous downloads still require the
processed files to be publicly released. Filenames are resolved from metadata,
not guessed from the accession.''')
s = replace_once(s,
    'python scripts/prepare_geo_plot_data.py --dataset hippo_timecourse --geo-accession GSE123456\npython scripts/prepare_geo_plot_data.py --dataset ipsc_k562_mix --local-file data/geo/iPSC_K562_mixed_integrated.h5mu',
    'python scripts/prepare_geo_plot_data.py --dataset hippo_timecourse --geo-accession GSE352939\npython scripts/prepare_geo_plot_data.py --dataset ipsc_k562_mix --geo-accession GSE352893')
s = replace_once(s,
    'without institutional paths. Actual accession-specific downloads require the\npublic study accession; they cannot be certified before it is available.',
    'without institutional paths. Accession-specific end-to-end regeneration requires\npublic access to the submitted processed objects; a configured accession or a\npassing mocked-download test does not establish that GEO has released those files.')
path.write_text(s)

Path('data/geo/README.md').write_text('''# Optional processed GEO objects

Default figure execution uses the compact tables in
`data/reproduction/geo_plots/`; no GEO download is required.

| Dataset | GEO accession | Notebook / panels |
|---|---|---|
| Parse v3 CC-Perturb-seq, iPSC-K562 mix | [GSE352893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352893) | Figure S2i, j, l, m, n |
| Small molecule Hippo scRNA-Seq | [GSE352939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352939) | Figures 3f and 4c |

These accessions are prefilled in the corresponding notebooks. Set
`DATA_SOURCE="geo"` to resolve and download the processed supplementary file and
regenerate plotting tables from the submitted integrated H5MU. `GEO_FILENAME`
may be set to an exact filename if discovery reports ambiguous candidates.

Downloads are cached in `data/geo/downloads/` and remain ignored by Git.
Regenerated tables go to `data/derived/geo_plots/`, never over published tables.
Released GEO files need no login. Before release, select `local_file` and provide
`Hippo_timecourse_integrated.h5mu` or `iPSC_K562_mixed_integrated.h5mu` from the
submission workflow. Assigned accessions do not by themselves guarantee public
file availability. Regeneration does not download FASTQs, refit cNMF, or
recluster cells.
''')

path = Path('tests/test_geo_plot_notebooks.py')
s = path.read_text()
s = replace_once(s,
    '.replace(\'GEO_ACCESSION = ""\',\'GEO_ACCESSION = "GSE12345"\')', '')
s = s.replace('GSE12nnn/GSE12345/suppl/GSE12345_iPSC_K562_mixed_integrated.h5mu',
              'GSE352nnn/GSE352893/suppl/GSE352893_iPSC_K562_mixed_integrated.h5mu')
s = replace_once(s, '^SERIES = GSE12345', '^SERIES = GSE352893')
s = replace_once(s,
    "    if 'acc.cgi?' in requested:\n        return Response(",
    "    if 'acc.cgi?' in requested:\n        from urllib.parse import urlparse,parse_qs\n        assert parse_qs(urlparse(requested).query)['acc']==['GSE352893']\n        return Response(")
s = replace_once(s,
    '''        assert 'GEO_ACCESSION = ""' in text''',
    '''        expected = 'GSE352893' if path.name == 'Figure_S02.ipynb' else 'GSE352939'
        assert f'GEO_ACCESSION = "{expected}"' in text''')
path.write_text(s)

Path('/tmp/geo-accession-preservation.json').write_text(json.dumps(reports, indent=2)+'\n')
print(json.dumps(reports, indent=2))

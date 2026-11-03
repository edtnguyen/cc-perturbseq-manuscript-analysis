"""Study accession wiring tests; mocked metadata is not a live GEO download."""
import ast
import io
import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest
from cc_perturbseq.geo_download import GeoClient

ROOT = Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('path,accession', [
    ('notebooks/main/Figure_03.ipynb', 'GSE352939'),
    ('notebooks/main/Figure_04.ipynb', 'GSE352939'),
    ('notebooks/supplementary/Figure_S02.ipynb', 'GSE352893'),
])
def test_notebook_option_dictionary_contains_correct_accession(path, accession):
    nb = json.loads((ROOT / path).read_text())
    variables = {'DATA_SOURCE', 'GEO_ACCESSION', 'GEO_FILENAME',
                 'GEO_EXPECTED_SHA256', 'GEO_LOCAL_FILE', 'GEO_OPTIONS'}
    statements = []
    for cell in nb['cells']:
        if cell['cell_type'] != 'code':
            continue
        for node in ast.parse(''.join(cell['source'])).body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in variables for t in node.targets):
                statements.append(node)
    namespace = {}
    exec(compile(ast.Module(body=statements, type_ignores=[]), path, 'exec'), namespace)
    assert namespace['DATA_SOURCE'] == 'precomputed'
    assert namespace['GEO_ACCESSION'] == accession
    assert namespace['GEO_OPTIONS'] == {
        'mode': 'precomputed', 'accession': '', 'filename': None,
        'expected_sha256': None, 'local_file': None,
    }
    # Changing the mode alone must retain the correct accession in GEO_OPTIONS.
    for node in statements:
        if any(isinstance(t, ast.Name) and t.id == 'DATA_SOURCE' for t in node.targets):
            node.value = ast.copy_location(ast.Constant(value='geo'), node.value)
    namespace = {}
    exec(compile(ast.Module(body=statements, type_ignores=[]), path, 'exec'), namespace)
    assert namespace['GEO_OPTIONS']['mode'] == 'geo'
    assert namespace['GEO_OPTIONS']['accession'] == accession

@pytest.mark.parametrize('dataset,accession,name', [
    ('ipsc_k562_mix', 'GSE352893', 'iPSC_K562_mixed_integrated.h5mu'),
    ('hippo_timecourse', 'GSE352939', 'Hippo_timecourse_integrated.h5mu'),
])
def test_study_accession_is_used_in_metadata_request(dataset, accession, name):
    # Synthetic attachment metadata, not a claim about current deposited names.
    supplementary = f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE352nnn/{accession}/suppl/{accession}_{name}'
    calls = []
    class FixtureClient(GeoClient):
        def _open(self, url):
            calls.append(url)
            query = parse_qs(urlparse(url).query)
            assert query['acc'] == [accession]
            return io.BytesIO((f'^SERIES = {accession}\n!Series_supplementary_file = {supplementary}\n').encode())
    record = FixtureClient().discover(accession, dataset=dataset)
    assert record['accession'] == accession
    assert record['dataset'] == dataset
    assert record['url'] == supplementary
    assert len(calls) == 1


def test_public_readmes_document_both_accessions_and_no_placeholders():
    for name in ('README.md', 'docs/PUBLIC_DATA.md', 'data/geo/README.md'):
        text = (ROOT / name).read_text()
        assert 'GSE352893' in text and 'GSE352939' in text
        assert 'GSE123456' not in text
        assert 'download URLs are not configured yet' not in text
        assert 'Public study accessions have not yet been supplied' not in text

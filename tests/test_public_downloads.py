from pathlib import Path
import hashlib
import io
import json
import os
from urllib.error import HTTPError
from urllib.request import Request
import pytest
from cc_perturbseq.igvf_download import API, IgvfClient, SafeRedirect
from cc_perturbseq.public_data import load_public_paths, require_bundle


class Opener:
    def __init__(self, record, payload, auth=False):
        self.record, self.payload, self.auth, self.requests = record, payload, auth, []
    def open(self, request, timeout=None):
        self.requests.append(request)
        if self.auth and not request.get_header('Authorization'):
            raise HTTPError(request.full_url,401,'denied',{},None)
        if request.full_url.endswith('?format=json'):
            return io.BytesIO(json.dumps({**self.record,'href':'/matrix-files/IGVFFI0000TEST/@@download/IGVFFI0000TEST.h5ad','upload_status':'validated'}).encode())
        return io.BytesIO(self.payload)


def record(payload=b'correct input'):
    return {'accession':'IGVFFI0000TEST','filename':'IGVFFI0000TEST.h5ad','md5sum':hashlib.md5(payload).hexdigest(),'file_size':len(payload)}


@pytest.fixture(autouse=True)
def no_real_credentials(monkeypatch):
    monkeypatch.delenv('IGVF_ACCESS_KEY',raising=False)
    monkeypatch.delenv('IGVF_SECRET_ACCESS_KEY',raising=False)


def test_public_download_checksum_and_cached_reuse(tmp_path):
    rec=record();op=Opener(rec,b'correct input');c=IgvfClient(allow_prompt=False,opener=op)
    p=c.download(rec,tmp_path)
    assert p.read_bytes()==b'correct input' and len(op.requests)==2
    assert c.download(rec,tmp_path)==p and len(op.requests)==2
    assert all(r.get_header('Authorization') is None for r in op.requests)
    assert not list(tmp_path.glob('*.part'))


def test_corrupt_download_does_not_replace_existing_cache(tmp_path):
    rec=record();p=tmp_path/rec['filename'];p.write_bytes(b'old')
    c=IgvfClient(allow_prompt=False,opener=Opener(rec,b'X'*rec['file_size']))
    with pytest.raises(ValueError,match='MD5'):
        c.download(rec,tmp_path)
    assert p.read_bytes()==b'old' and not list(tmp_path.glob('*.part'))


def test_pin_mismatch_stops_before_downloading(tmp_path):
    rec=record();op=Opener({**rec,'md5sum':'0'*32},b'correct input')
    with pytest.raises(ValueError,match='pinned md5sum'):
        IgvfClient(allow_prompt=False,opener=op).download(rec,tmp_path)
    assert len(op.requests)==1 and not list(tmp_path.glob('*.h5ad'))


def test_auth_prompt_is_hidden_optional_and_not_persisted(tmp_path,monkeypatch,capsys):
    answers=iter(['TESTKEY','TESTSECRET']);monkeypatch.setattr('getpass.getpass',lambda _:next(answers))
    op=Opener(record(),b'correct input',auth=True);c=IgvfClient(opener=op)
    c.download(record(),tmp_path)
    assert any(r.get_header('Authorization') for r in op.requests)
    assert 'TESTSECRET' not in repr(c)+capsys.readouterr().out
    for p in tmp_path.glob('*.json'):
        assert 'TESTKEY' not in p.read_text() and 'TESTSECRET' not in p.read_text()


def test_unauthorized_batch_never_prompts(tmp_path,monkeypatch):
    monkeypatch.setattr('getpass.getpass',lambda _:pytest.fail('unexpected prompt'))
    with pytest.raises(PermissionError):
        IgvfClient(allow_prompt=False,opener=Opener(record(),b'',auth=True)).download(record(),tmp_path)


def test_redirect_drops_auth_on_storage_origin_and_rejects_http():
    req=Request(API+'/x',headers={'Authorization':'Basic TESTSECRET'})
    out=SafeRedirect().redirect_request(req,None,302,'redirect',{},'https://bucket.s3.amazonaws.com/data?signature=test')
    assert out.get_header('Authorization') is None
    same=SafeRedirect().redirect_request(req,None,302,'redirect',{},API+'/next')
    assert same.get_header('Authorization')=='Basic TESTSECRET'
    with pytest.raises(RuntimeError,match='non-HTTPS'):
        SafeRedirect().redirect_request(req,None,302,'redirect',{},'http://bucket.test/x')


def test_manifest_cannot_direct_authenticated_request_to_other_origin(tmp_path):
    with pytest.raises(ValueError,match='official HTTPS'):
        IgvfClient()._open('https://attacker.test/file')


def test_legacy_local_config_does_not_override_portable_paths(tmp_path,monkeypatch):
    (tmp_path/'config').mkdir();(tmp_path/'config/paths.example.yaml').write_text('outputs:\n  root: output\n')
    (tmp_path/'config/paths.local.yaml').write_text('outputs:\n  root: /oak/private/path\n')
    monkeypatch.delenv('CCPERTURBSEQ_PUBLIC_PATHS',raising=False)
    with pytest.warns(UserWarning,match='Legacy'):
        _,paths=load_public_paths(tmp_path)
    assert paths['outputs']['root']==str(tmp_path/'output')


def test_explicit_institutional_override_is_rejected(tmp_path,monkeypatch):
    (tmp_path/'config').mkdir();(tmp_path/'config/paths.example.yaml').write_text('outputs:\n  root: output\n')
    override=tmp_path/'override.yaml';override.write_text('outputs:\n  root: /oak/somewhere/else\n')
    monkeypatch.setenv('CCPERTURBSEQ_PUBLIC_PATHS',str(override))
    with pytest.raises(ValueError,match='institutional'):
        load_public_paths(tmp_path)


def test_companion_checksum_is_enforced(tmp_path):
    data=tmp_path/'data/reproduction';data.mkdir(parents=True)
    p=data/'input.tsv';p.write_text('scientific inputs\n')
    manifest=data/'bundle_manifest.json';manifest.write_text(json.dumps({'files':[{'key':'test','path':'data/reproduction/input.tsv','size_bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}]}))
    paths={'public_data':{'bundle_manifest':str(manifest)}}
    assert require_bundle(tmp_path,paths,['test'])['test']==p
    p.write_text('different values\n')
    with pytest.raises(ValueError,match='checksum'):
        require_bundle(tmp_path,paths,['test'])


def test_missing_bundle_is_not_faked(tmp_path):
    with pytest.raises(FileNotFoundError,match='not been published'):
        require_bundle(tmp_path,{'public_data':{'bundle_manifest':str(tmp_path/'missing.json')}},['test'])

import gzip
import hashlib
import io
import json
from pathlib import Path
import h5py
import pytest
from cc_perturbseq.geo_download import GeoClient, ncbi_url, normalize_accession

URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12345/suppl/GSE12345_Hippo_timecourse_integrated.h5mu'

class Response(io.BytesIO):
    def __init__(self, raw, length=None):
        super().__init__(raw)
        self.headers = {'Content-Length': str(len(raw) if length is None else length)}

class Client(GeoClient):
    def __init__(self, records, files):
        self.records = records; self.files = files; self.calls=[]
    def _open(self, url):
        from urllib.parse import urlparse, parse_qs
        self.calls.append(url)
        query=parse_qs(urlparse(url).query)
        if 'acc' in query:
            return Response(self.records[query['targ'][0]].encode())
        value=self.files[url]
        return value() if callable(value) else Response(value)


def metadata(*urls, accession='GSE12345', kind='Series'):
    return '^'+('SERIES' if kind=='Series' else 'SAMPLE')+' = '+accession+'\n'+''.join('!'+kind+'_supplementary_file = '+url+'\n' for url in urls)


def h5bytes(tmp_path):
    path=tmp_path/'valid.h5mu'
    with h5py.File(path,'w') as f:f.create_group('mod')
    return path.read_bytes()

@pytest.mark.parametrize('acc',['','SRA1234','../GSE1234','GSE1234; rm','https://example.com','GSM0'])
def test_invalid_accession_rejected(acc):
    with pytest.raises(ValueError):normalize_accession(acc)


def test_gse_and_gsm_metadata_discovery_and_ftp_upgrade():
    assert normalize_accession(' gse12345 ')=='GSE12345'
    client=Client({'self':metadata(URL.replace('https:','ftp:'))},{})
    assert client.discover('GSE12345',dataset='hippo_timecourse')['url']==URL
    client=Client({'self':metadata(URL,accession='GSM12345',kind='Sample')},{})
    assert client.discover('GSM12345',dataset='hippo_timecourse')['url']==URL


def test_sample_level_discovery_and_ambiguous_selection():
    second=URL.replace('integrated','other')
    client=Client({'self':metadata(), 'gsm':metadata(URL,second,accession='GSM9999',kind='Sample')},{})
    with pytest.raises(ValueError,match='found 2'):client.discover('GSE12345',dataset='hippo_timecourse')
    assert client.discover('GSE12345',dataset='hippo_timecourse',filename=Path(URL).name)['url']==URL


def test_tar_is_not_substituted_and_nonpublic_metadata_fails():
    client=Client({'self':metadata(URL+'.tar'),'gsm':''},{})
    with pytest.raises(ValueError,match='found 0'):client.discover('GSE12345',dataset='hippo_timecourse')
    client.records['self']='<html>record is not public</html>'
    with pytest.raises(FileNotFoundError,match='public SOFT'):client.discover('GSE12345',dataset='hippo_timecourse')

@pytest.mark.parametrize('url',['https://evil.example/file.h5mu','file:///etc/passwd','http://ftp.ncbi.nlm.nih.gov/geo/a','https://user:secret@ftp.ncbi.nlm.nih.gov/geo/a','https://ftp.ncbi.nlm.nih.gov/not-geo/a'])
def test_unsafe_hosts_schemes_or_credentials_rejected(url):
    with pytest.raises(ValueError):ncbi_url(url)


def test_real_hdf5_payload_gzip_checksum_cache_and_partial_safety(tmp_path):
    raw=h5bytes(tmp_path);wire=gzip.compress(raw);url=URL+'.gz'
    client=Client({'self':metadata(url)}, {url:wire})
    expected=hashlib.sha256(raw).hexdigest()
    path=client.download('GSE12345',tmp_path/'cache',dataset='hippo_timecourse',expected_sha256=expected)
    assert path.read_bytes()==raw
    assert json.loads(path.with_name(path.name+'.download.json').read_text())['expected_sha256_verified']
    client.download('GSE12345',tmp_path/'cache',dataset='hippo_timecourse',expected_sha256=expected)
    assert client.calls.count(url)==1
    # Corruption is replaced, not reused as a valid cache hit.
    path.write_bytes(b'bad')
    client.download('GSE12345',tmp_path/'cache',dataset='hippo_timecourse',expected_sha256=expected)
    assert path.read_bytes()==raw and client.calls.count(url)==2
    client.files[url]=lambda:Response(wire,length=len(wire)+5)
    with pytest.raises(ValueError,match='Incomplete'):
        client.download('GSE12345',tmp_path/'partial',dataset='hippo_timecourse')
    assert not list((tmp_path/'partial').rglob('*.h5mu'))


def test_wrong_hash_and_html_never_published(tmp_path):
    raw=h5bytes(tmp_path);client=Client({'self':metadata(URL)}, {URL:raw})
    with pytest.raises(ValueError,match='SHA-256 mismatch'):
        client.download('GSE12345',tmp_path/'cache',dataset='hippo_timecourse',expected_sha256='0'*64)
    client.files[URL]=b'<html>not a matrix</html>'
    with pytest.raises(ValueError,match='container'):
        client.download('GSE12345',tmp_path/'cache',dataset='hippo_timecourse')
    assert not list((tmp_path/'cache').rglob('*.h5mu'))

from pathlib import Path
import base64, hashlib, json, lzma, subprocess, tempfile
BASE='8c70476231298d41e08184e75b687bddb524a52a'
SHA='b941e94a9a251603066ef669c58ed9550ea22dfe3f9b92a15f022d4a62776ed3'
def digest(x):return hashlib.sha256(x).hexdigest()
def preserved(nb):
    n=json.loads(json.dumps(nb))
    for c in n['cells']:c['source']=None
    return digest(json.dumps(n,sort_keys=True,ensure_ascii=False).encode())
def patch(raw,text):
    with tempfile.TemporaryDirectory() as folder:
        p=Path(folder);(p/'cell').write_bytes(raw);(p/'patch').write_text(text)
        subprocess.run(['patch','--batch','--fuzz=0','-p1','-i','patch'],cwd=p,check=True)
        return (p/'cell').read_bytes()
def main():
    import sys
    raw=lzma.decompress(base64.b64decode(Path(sys.argv[1]).read_text().strip()))
    assert digest(raw)==SHA,'Payload checksum mismatch'
    payload=json.loads(raw)
    assert payload['baseline']==BASE
    assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==BASE
    paths=[]
    for spec in payload['files']:
        p=Path(spec['path']);assert not p.is_absolute() and '..' not in p.parts
        if 'patch' in spec:
            old=p.read_bytes();assert digest(old)==spec['old_sha256'],str(p)
            new=patch(old,spec['patch']);assert digest(new)==spec['new_sha256'],str(p)
        else:
            assert not p.exists(),str(p);new=spec['new'].encode()
        p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(new);paths.append(str(p))
    reports=[]
    for spec in payload['notebooks']:
        p=Path(spec['path']);nb=json.loads(p.read_text());before=preserved(nb)
        assert len(nb['cells'])==spec['cell_count']
        for change in spec['cells']:
            c=nb['cells'][change['index']];old=''.join(c['source']).encode()
            assert digest(old)==change['old_sha256'],(p,change['index'])
            new=patch(old,change['patch']);assert digest(new)==change['new_sha256']
            c['source']=new.decode().splitlines(keepends=True)
        assert preserved(nb)==before,str(p)
        for i,c in enumerate(nb['cells']):
            if c['cell_type']=='code':compile(''.join(c['source']),f'{p}:{i}','exec')
        p.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')
        assert preserved(json.loads(p.read_text()))==before,str(p)
        reports.append({'path':str(p),'unchanged_non_source_sha256':before,'outputs_preserved':sum(len(c.get('outputs',[])) for c in nb['cells'])})
        paths.append(str(p))
    Path('/tmp/output-env-paths.json').write_text(json.dumps(paths))
    Path('/tmp/output-env-preservation.json').write_text(json.dumps(reports,indent=2))
    for p in Path('src').rglob('*.py'):compile(p.read_text(),str(p),'exec')
    print('Applied',len(paths),'files; preserved every notebook output and metadata field.')
if __name__=='__main__':main()

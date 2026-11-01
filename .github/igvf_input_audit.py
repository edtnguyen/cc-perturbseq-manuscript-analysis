"""Read-only public-data audit. No portal credentials and no changes to main."""
import concurrent.futures
import datetime
import gzip
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

OUT=Path('/tmp/igvf-input-audit')
OUT.mkdir(parents=True,exist_ok=True)
BASE='https://api.data.igvf.org'
SETS=['IGVFDS9629VXQK','IGVFDS4708PCVX','IGVFDS0463XTUH','IGVFDS0991DBQZ','IGVFDS4402ZGWJ','IGVFDS2307REVJ','IGVFDS5257JGQO','IGVFDS3530HHEU','IGVFDS7380DUFB','IGVFDS6826YKCZ','IGVFDS9440URJB','IGVFDS9710KKQV']
EXTRA_FILES=['IGVFFI0180ZBTX','IGVFFI3590IELV','IGVFFI6391GCJR','IGVFFI6745HXWN','IGVFFI9629WVUU']

def request(url,method='GET'):
    return urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'Manuscript-input-lineage-audit'},method=method)

def get(path):
    url=BASE+'/'+path.lstrip('/')
    try:
        with urllib.request.urlopen(request(url),timeout=50) as r:
            return {'url':url,'status_code':r.status,'data':json.load(r)}
    except Exception as exc:
        return {'url':url,'error':str(exc)}

def save(key,value):
    (OUT/(key+'.json')).write_text(json.dumps(value,indent=2,default=str)+'\n')

def ref(x):
    return x.get('accession',x.get('@id')) if isinstance(x,dict) else x

records={}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for key,item in zip(SETS,pool.map(lambda a:get(a+'/?format=json'),SETS)):
        records[key]=item;save(key,item)
files=set(EXTRA_FILES)
for item in records.values():
    for f in item.get('data',{}).get('files',[]):
        value=ref(f)
        if value:files.update(re.findall(r'IGVFFI[A-Z0-9]{8}',value))
for depth in range(3):
    missing=sorted(files-set(records))
    if not missing:break
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for key,item in zip(missing,pool.map(lambda a:get(a+'/?format=json'),missing)):
            records[key]=item;save(key,item)
            for f in item.get('data',{}).get('derived_from',[]):
                value=ref(f)
                if value and '/sequence-files/' not in value and '/alignment-files/' not in value:
                    files.update(re.findall(r'IGVFFI[A-Z0-9]{8}',value))

manifest=[]
for accession,record in sorted(records.items()):
    if not accession.startswith('IGVFFI'):continue
    d=record.get('data',{})
    entry={k:d.get(k) for k in ['accession','@id','content_type','file_format','status','upload_status','submitted_file_name','md5sum','file_size','href','s3_uri','description','derived_from','aliases']}
    entry['file_set']=ref(d.get('file_set'))
    manifest.append(entry)
save('file_manifest',manifest)

# Download public tabular files only, with a per-file limit. Matrices are
# inspected separately; only their schema/annotations are retained as artifacts.
DOWNLOADS=OUT/'tables';DOWNLOADS.mkdir(exist_ok=True)
def download_table(entry):
    acc=entry['accession']
    if entry.get('file_format') not in {'tsv','csv'} or not entry.get('href'):
        return None
    if not entry.get('file_size') or entry['file_size']>40_000_000:
        return {'accession':acc,'skipped':'size limit'}
    import pandas as pd
    try:
        url=BASE+entry['href']
        filename=Path(urllib.parse.urlparse(entry['href']).path).name
        path=DOWNLOADS/filename
        with urllib.request.urlopen(request(url),timeout=90) as response,path.open('wb') as target:
            while True:
                block=response.read(1024*1024)
                if not block:break
                target.write(block)
        observed=hashlib.md5(path.read_bytes()).hexdigest()
        if observed!=entry['md5sum']:raise ValueError('MD5 mismatch')
        frame=pd.read_csv(path,sep='\t' if entry['file_format']=='tsv' else ',',low_memory=False)
        return {'accession':acc,'path':str(path.relative_to(OUT)),'md5_matches':True,
                'rows':len(frame),'columns':list(frame.columns),'preview':frame.head(3).fillna('').to_dict('records')}
    except Exception as exc:
        return {'accession':acc,'error':str(exc)}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    summaries=[v for v in pool.map(download_table,manifest) if v]
save('table_schemas',summaries)

# h5py over HTTP Range requests avoids loading count matrices. If the public
# server does not support Range, a bounded local download is used instead.
def inspect_matrix(accession):
    import h5py
    import numpy as np
    import fsspec
    d=records[accession]['data']
    url=BASE+d['href']
    report={'accession':accession,'submitted_file_name':d['submitted_file_name'],'file_size':d['file_size']}
    stream=None;handle=None
    try:
        with urllib.request.urlopen(request(url,method='HEAD'),timeout=60) as response:
            resolved=response.url
        try:
            stream=fsspec.open(resolved,mode='rb',block_size=1024*1024).open()
            handle=h5py.File(stream,'r')
            report['access_method']='HTTP Range'
        except Exception as exc:
            if stream:stream.close()
            if d['file_size']>1_000_000_000:raise
            local=Path('/tmp')/(accession+'.'+d['file_format'])
            with urllib.request.urlopen(request(url),timeout=120) as response,local.open('wb') as target:
                while True:
                    block=response.read(4*1024*1024)
                    if not block:break
                    target.write(block)
            md5=hashlib.md5()
            with local.open('rb') as f:
                for block in iter(lambda:f.read(4*1024*1024),b''):md5.update(block)
            if md5.hexdigest()!=d['md5sum']:raise ValueError('Matrix download MD5 mismatch')
            handle=h5py.File(local,'r');report['access_method']='local download with matching MD5'
        datasets={};groups={}
        def convert(x):
            if isinstance(x,bytes):return x.decode(errors='replace')
            if isinstance(x,np.ndarray):return [convert(v) for v in x.tolist()]
            if isinstance(x,np.generic):return convert(x.item())
            if isinstance(x,list):return [convert(v) for v in x]
            return x
        def visit(name,obj):
            if isinstance(obj,h5py.Group):
                groups[name]={'keys':list(obj.keys()),'attrs':{k:convert(v) for k,v in obj.attrs.items()}}
                return
            value={'shape':list(obj.shape),'dtype':str(obj.dtype),'compression':obj.compression}
            try:
                if obj.size<=1000:
                    value['value']=convert(obj[()])
                elif obj.ndim==1 and ('/var/' in '/'+name or '/obs/' in '/'+name or name.startswith('uns/')):
                    value['preview']=convert(obj[:10])
                elif obj.ndim==2 and obj.dtype.kind in 'fiu':
                    value['preview']=convert(obj[:3,:min(5,obj.shape[1])])
            except Exception as exc:value['read_error']=str(exc)
            datasets[name]=value
        handle.visititems(visit)
        report['groups']=groups;report['datasets']=datasets
        # Preserve lightweight full program annotations, gene IDs, and usage row
        # totals for identifying semantics; do not export single-cell matrices.
        for prefix in ['mod/cNMF','mod/rna','']:
            if prefix and prefix not in handle:continue
            group=handle[prefix] if prefix else handle
            if 'var' in group:
                vg=group['var'];index_name=vg.attrs.get('_index','_index')
                if isinstance(index_name,bytes):index_name=index_name.decode()
                if index_name in vg:
                    report.setdefault('var_indices',{})[prefix or 'root']=convert(vg[index_name][()])
            if prefix=='mod/cNMF' and 'X' in group and isinstance(group['X'],h5py.Dataset):
                x=group['X'][:10,:]
                report['cnmf_first_10_row_sums']=convert(np.sum(x,axis=1))
        report['status']='inspected'
    except Exception as exc:
        report.update(status='failed',error=str(exc))
    finally:
        if handle:handle.close()
        if stream:stream.close()
    save(accession+'_hdf5_schema',report)
    print(accession,report['status'],report.get('error',''))

for accession in ['IGVFFI8722NAWN','IGVFFI7756RKHT']:
    inspect_matrix(accession)
save('audit_context',{'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'analysis_sets':SETS,'scope':'Read-only public metadata, TSV payloads and two representative HDF5 schemas; main unchanged'})
print('Audited',len(manifest),'files;',len(summaries),'table downloads')

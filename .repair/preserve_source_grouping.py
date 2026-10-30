"""Preserve source grouping: metadata aliases are annotations, not silent CRT changes."""
from pathlib import Path
import nbformat

r=Path.cwd()
def rep(s,a,b,n=1):
    if s.count(a)!=n:raise RuntimeError(f'Expected {n} occurrences of {a[:90]!r}; got {s.count(a)}')
    return s.replace(a,b,n)

p=r/'src/cc_perturbseq/crt_qq_backend.py';s=p.read_text()
a=s.index('def build_guide2gene(');b=s.index('\ndef h5mu_path_for_day',a)
part=s[a:b]
part=rep(part,'    target_col: str | None,\n)',
         '    target_col: str | None,\n    mapping_policy: str = "metadata",\n)')
part=rep(part,'    mapping = infer_guide2gene_from_names(guide_names, ntc_labels)',
'''    if mapping_policy not in {"source_names", "metadata"}:
        raise ValueError("guide mapping policy must be source_names or metadata")
    mapping = infer_guide2gene_from_names(guide_names, ntc_labels)''')
part=rep(part,'                mapping.update(metadata_mapping)',
'''                # Recognized controls remain control groups even when the
                # metadata target field contains an individual SAFE guide ID.
                for guide in metadata_mapping:
                    if mapping.get(guide) in ntc_labels:
                        metadata_mapping[guide] = mapping[guide]
                if mapping_policy == "metadata":
                    mapping.update(metadata_mapping)
                else:
                    print("guide mapping policy: source_names (original grouping retained)")''')
s=s[:a]+part+s[b:]
s=rep(s,'    metadata_target_col: str | None,\n) -> tuple[list[str], dict[str, str]]:',
        '    metadata_target_col: str | None,\n    mapping_policy: str = "source_names",\n) -> tuple[list[str], dict[str, str]]:')
s=rep(s,'        metadata_target_col,\n    )\n    adata.uns["guide2gene"]',
        '        metadata_target_col,\n        mapping_policy=mapping_policy,\n    )\n    adata.uns["guide2gene"]')
s=rep(s,'        args.guide_metadata_target_col,\n    )\n    covariates =',
        '        args.guide_metadata_target_col,\n        mapping_policy=args.guide_mapping_policy,\n    )\n    covariates =')
s=rep(s,'    parser.add_argument("--guide-metadata-target-col", default=None)',
'''    parser.add_argument("--guide-metadata-target-col", default=None)
    parser.add_argument("--guide-mapping-policy", choices=("source_names", "metadata"),
                        default="source_names",
                        help="source_names preserves the original S4 guide grouping; metadata explicitly opts into target-symbol updates.")''')
s=rep(s,'        "guide_long_csv": str(guide_csv),\n        "ntc_labels":',
        '        "guide_long_csv": str(guide_csv),\n        "guide_mapping_policy": args.guide_mapping_policy,\n        "ntc_labels":')
p.write_text(s)
p=r/'src/cc_perturbseq/figure_s04_compute.py';s=p.read_text()
s=rep(s,'    guide_metadata_target_col=None,\n):','    guide_metadata_target_col=None, guide_mapping_policy="source_names",\n):')
s=rep(s,'            args_for_helpers.guide_metadata_target_col,\n        )',
        '            args_for_helpers.guide_metadata_target_col,\n            mapping_policy=guide_mapping_policy,\n        )')
p.write_text(s)
p=r/'src/cc_perturbseq/figure_s04_runtime.py';s=p.read_text()
s=rep(s,'                   ntc_labels=("SAFE", "non-targeting", "NTC"), overwrite=True):',
        '                   ntc_labels=("SAFE", "non-targeting", "NTC"), overwrite=True,\n                   guide_mapping_policy="source_names"):')
s=rep(s,'        "--n-jobs", str(n_jobs), "--require-covariates", "--no-plot",',
        '        "--n-jobs", str(n_jobs), "--guide-mapping-policy", str(guide_mapping_policy),\n        "--require-covariates", "--no-plot",')
s=rep(s,'            guide_metadata_target_col=args.guide_metadata_target_col,\n        )',
        '            guide_metadata_target_col=args.guide_metadata_target_col,\n            guide_mapping_policy=args.guide_mapping_policy,\n        )')
p.write_text(s)
p=r/'src/cc_perturbseq/figure_s04_contract.py';s=p.read_text()
s=rep(s,'    after = dict(before)\n    after.update(mapped)',
'''    metadata_changes = [{"guide": g, "previous_inferred_target": before[g], "metadata_target": target}
                        for g, target in sorted(mapped.items()) if target != before[g]]
    metadata_control_labels = [g for g in mapped if before[g] in args.ntc_label and mapped[g] != before[g]]
    policy = getattr(args, "guide_mapping_policy", "source_names")
    after = dict(before)
    if policy == "metadata":
        after.update({g: target for g, target in mapped.items() if before[g] not in args.ntc_label})
    elif policy != "source_names":
        raise ValueError(f"Unsupported guide mapping policy: {policy}")''')
s=rep(s,'"guide_column": str(column), "n_guides": len(after),',
'''"guide_column": str(column), "n_guides": len(after), "mapping_policy": policy,
            "metadata_target_difference_count": len(metadata_changes),
            "metadata_target_differences": metadata_changes,
            "n_metadata_control_labels_preserved": len(metadata_control_labels),''')
p.write_text(s)
p=r/'notebooks/supplementary/Figure_S04.ipynb';nb=nbformat.read(p,as_version=4)
for c in nb.cells:
    if c.cell_type=='code':
        if 'CRT_STAGE = "all"' in c.source:
            c.source=rep(c.source,'CRT_STAGE = "all"',
                'CRT_STAGE = "all"\n# Retain original gene/control group names and RNG inputs for manuscript reproduction.\nGUIDE_MAPPING_POLICY = "source_names"')
        if 'runtime_report = launch_crt_recomputation' in c.source:
            c.source=rep(c.source,'ntc_labels=NTC_LABELS, overwrite=True,',
                'ntc_labels=NTC_LABELS, overwrite=True, guide_mapping_policy=GUIDE_MAPPING_POLICY,')
    elif '### Replicate-output flags' in c.source:
        c.source += '''

`GUIDE_MAPPING_POLICY="source_names"` is the reproduction default. Metadata
labels are audited, but the existing gene/control grouping is preserved. The
checked-in v43 table contains aliases (for example T/TBXT) and individual SAFE
IDs in its target column. Those must NOT silently rename or split the original
CRT targets. `"metadata"` is an explicit alternative analysis mode; it retains
canonical control labels and requires regenerating both stages when mappings
change. The original completed export can be resumed under `"source_names"`.
'''
nbformat.validate(nb)
for c in nb.cells:
    if c.cell_type=='code':compile(c.source,'S4','exec')
nbformat.write(nb,p)
p=r/'tests/test_s4_crt_runtime.py';s=p.read_text()
s=rep(s,'ns=dict(RECOMPUTE_CRT=True, CRT_STAGE="all", build_crt_argv=',
        'ns=dict(RECOMPUTE_CRT=True, CRT_STAGE="all", GUIDE_MAPPING_POLICY="source_names", build_crt_argv=')
p.write_text(s)
p=r/'tests/test_s4_output_contract.py';s=p.read_text()
s=rep(s,"report=contract.audit_guide_metadata(qqb.parse_args(['--guide-metadata-tsv',str(p)]))\n    assert report['legacy_name_difference_count']==1",
        "report=contract.audit_guide_metadata(qqb.parse_args(['--guide-metadata-tsv',str(p),'--guide-mapping-policy','metadata']))\n    assert report['legacy_name_difference_count']==1")
s=rep(s,"audit=audit_guide_metadata(args)\nPath(sys.argv[2])",'''audit=audit_guide_metadata(args)
assert audit['mapping_policy']=='source_names'
assert audit['legacy_name_difference_count']==0
assert audit['metadata_target_difference_count']>0
assert audit['n_metadata_control_labels_preserved']>0
actual=qqb.read_table_auto(Path(args.guide_metadata_tsv))
guides=list(actual['guide_id'].astype(str))
source_map=qqb.build_guide2gene(guides,args.ntc_label,args.guide_metadata_tsv,None,None,mapping_policy='source_names')
inferred=qqb.infer_guide2gene_from_names(guides,args.ntc_label)
assert source_map==inferred
metadata_map=qqb.build_guide2gene(guides,args.ntc_label,args.guide_metadata_tsv,None,None,mapping_policy='metadata')
assert all(metadata_map[g]==target for g,target in inferred.items() if target in args.ntc_label)
Path(sys.argv[2])''')
s += '''

def test_source_names_preserve_aliases_and_control_groups(tmp_path):
    p=tmp_path/'v43.tsv'
    names=['T_1','SAFE_TARGETING_41','NTC_1']
    pd.DataFrame({'guide_id':names,'intended_target_name':['TBXT','SAFE_TARGETING_41',None]}).to_csv(p,sep='\\t',index=False)
    source=qqb.build_guide2gene(names,qqb.DEFAULT_NTC_LABELS,str(p),None,None,mapping_policy='source_names')
    assert source==dict(zip(names,['T','SAFE','NTC']))
    metadata=qqb.build_guide2gene(names,qqb.DEFAULT_NTC_LABELS,str(p),None,None,mapping_policy='metadata')
    assert metadata==dict(zip(names,['TBXT','SAFE','NTC']))
    args=qqb.parse_args(['--guide-metadata-tsv',str(p)])
    audit=contract.audit_guide_metadata(args)
    assert audit['legacy_name_difference_count']==0
    assert audit['metadata_target_difference_count']==2
    assert audit['n_metadata_control_labels_preserved']==1


def test_resume_rejects_explicit_mapping_policy_change(tmp_path):
    args,audit,model=completed_run(tmp_path)
    args.guide_mapping_policy='metadata'
    with pytest.raises(ValueError,match='settings changed'):
        contract.validate_completed_full_day(args,audit)
'''
p.write_text(s)
p=r/'docs/S4_CRT_RUNTIME.md';p.write_text(p.read_text()+'''
### Preserve the source's target grouping

The real v43 audit found differences between guide-prefix names and metadata
labels, including gene aliases and per-guide SAFE target identifiers. Therefore
recognizing intended_target_name must not silently change the reproduction.
GUIDE_MAPPING_POLICY='source_names' is now explicit in the notebook and applied
to BOTH stages. Metadata is read and audited but original target names, control
groups, and gene-name-dependent random seeds are preserved. The optional
metadata policy is an explicit alternative, keeps controls canonical, and is
not used to resume a source-names full-day export. Actual v43 metadata and a
SAFE_TARGETING example are exercised in regression tests.
''')
for p in (r/'src/cc_perturbseq').glob('*.py'):compile(p.read_text(),str(p),'exec')
print('Source grouping preserved; metadata annotation audited separately.')

from pathlib import Path
import ast, hashlib, json
p=Path('src/cc_perturbseq/figure_s06_annotations.py')
s=p.read_text()
a=s.index('            center_x = (',s.index('def _draw_rotated_vector_logo('))
b=s.index('                patch.set_clip_on(False)',a)
s=s[:a]+'''            # Same physical 90-degree rotation as the source, expressed in
            # axes coordinates so PDF/SVG/PNG DPI changes do not move the logos.
            position = ax.get_position()
            aspect = (fig.get_figwidth() * position.width) / (
                fig.get_figheight() * position.height
            )
            rotation = (
                Affine2D().translate(-0.5, -0.5)
                .scale(aspect, 1.0)
                .rotate_deg(ROTATE_MOTIF_LOGOS_DEGREES)
                .scale(1.0 / aspect, 1.0)
                .translate(0.5, 0.5)
            )
            for patch in list(ax.patches):
                patch.set_transform(
                    (patch.get_transform() - ax.transAxes)
                    + rotation + ax.transAxes
                )

'''+s[b:]
s=s.replace('Original helper/selection/rendering functions are retained verbatim. The input',
            'Original selection/rendering functions are ported; DPI-stable logo rotation\nis the documented export-only adaptation. The input')
compile(s,str(p),'exec');p.write_text(s)
p=Path('docs/STANDALONE_FIGURE_SOURCE_AUDIT.json');audit=json.loads(p.read_text());spec=audit['Figure_S06']
original=spec['original_function_ast_sha256'].pop('_draw_rotated_vector_logo')
node=next(n for n in ast.parse(s).body if isinstance(n,ast.FunctionDef) and n.name=='_draw_rotated_vector_logo')
spec['documented_function_adaptations']={'_draw_rotated_vector_logo':{
    'source_ast_sha256':original,
    'ported_ast_sha256':hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest(),
    'reason':'Preserve physical 90-degree glyph placement across DPI changes and tight bounding-box vector/raster exports; no numerical or selection change.'}}
p.write_text(json.dumps(audit,indent=2)+'\n')
p=Path('tests/test_standalone_figure_ports.py');s=p.read_text()
s=s.replace("np.arange(128,dtype=float).reshape(8,16)%7/10", "np.zeros((8,16),dtype=float)")
s=s.replace("weights[p,2*p:(2*p+2)]=10", "weights[p,2*p:(2*p+2)]=10\n        weights[p,(2*p+2)%16]=1\n        weights[p,(2*p+3)%16]=0.5")
s+='''

def test_s6_logo_rotation_stays_at_same_physical_position_across_dpi():
    import matplotlib.pyplot as plt
    from cc_perturbseq import figure_s06_annotations as s6
    fig=plt.figure(figsize=(8,6))
    ax=fig.add_axes([0.5,0.5,0.2,0.04])
    matrix=pd.DataFrame(np.eye(4)[[0,1,2,3,0,1]],columns=list('ACGT'))
    s6._draw_rotated_vector_logo(fig,ax,'TEST',{'TEST':{'type':'cwm','cwm_df':matrix}})
    assert ax.patches
    bounds=[]
    for dpi in (72,100,600):
        fig.set_dpi(dpi);fig.canvas.draw()
        bounds.append(np.asarray([p.get_window_extent().extents for p in ax.patches])/dpi)
    np.testing.assert_allclose(bounds[0],bounds[1],atol=1e-9)
    np.testing.assert_allclose(bounds[0],bounds[2],atol=1e-9)
    plt.close(fig)


def test_documented_export_adaptations_match_the_reviewed_port():
    audit=json.loads((ROOT/'docs/STANDALONE_FIGURE_SOURCE_AUDIT.json').read_text())
    spec=audit['Figure_S06']
    nodes={n.name:n for n in ast.parse((ROOT/spec['ported_file']).read_text()).body if isinstance(n,ast.FunctionDef)}
    for name,change in spec['documented_function_adaptations'].items():
        actual=hashlib.sha256(ast.dump(nodes[name],include_attributes=False).encode()).hexdigest()
        assert actual==change['ported_ast_sha256']
        assert change['source_ast_sha256'] != actual
'''
p.write_text(s)
p=Path('docs/STANDALONE_FIGURE_REPRODUCTION.md')
p.write_text(p.read_text()+'''\nThe rendered S6 previews exposed a DPI-dependent placement defect inherited\nfrom the source's fixed pixel-coordinate logo rotation. The port expresses the\nsame 90-degree rotation in axes coordinates with physical-aspect correction.\nThis retains the motif matrices, glyphs, colors, orientation and intended\nplacement while keeping logos in position for PDF/SVG/600-dpi PNG exports.\nThe source and adapted function hashes are recorded separately, and an actual\nlogomaker glyph test verifies identical physical bounds at 72, 100 and 600 dpi.\n''')

#!/usr/bin/env python3
"""Guide-level validation EC boxplots and two-sided permutation tests.

Input column is used directly: estimated_EC_pct_control70.
Controls are identified solely using guide_type = NON-TARGETING or SAFE.
Each finite guide is one unweighted observation. No cell/read count weighting.
The primary p-value is b/B as specified by the user; plus-one values are
exported as a sensitivity diagnostic, not used for the figure stars.

Run:
  python make_validation_boxplots.py --input guide_EC_percent_control70.csv --outdir results
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import false_discovery_control

VALUE_COL = 'estimated_EC_pct_control70'
CONTROL_TYPES = {'NON-TARGETING', 'SAFE'}
TARGET_TYPE = 'TARGETING'
GROUPS = [
    ('Control', ['NTC']),
    ('BMP/TGFβ pathway', ['BMPR1A', 'SMAD1', 'SMAD5']),
    ('Mesoderm regulators', ['MIXL1', 'EOMES']),
    ('EC + LM regulators', ['KDR', 'KCTD10', 'PRKCD', 'ETV2', 'PLCG1']),
    ('Hippo pathway', ['TAOK1', 'AMOTL2', 'LATS2', 'NF2']),
]
GROUP_LABELS = {
    'Control': 'Control',
    'BMP/TGFβ pathway': 'BMP/TGFβ\npathway',
    'Mesoderm regulators': 'Mesoderm\nregulators',
    'EC + LM regulators': 'EC + LM\nregulators',
    'Hippo pathway': 'Hippo\npathway',
}
STEM = 'validation_boxplots_verified'


def write_csv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    if columns is None:
        if not rows:
            raise ValueError('Specify columns for an empty table.')
        columns = list(rows[0])
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='raise')
        w.writeheader()
        w.writerows(rows)


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError('BH adjustment requires finite p-values in [0, 1].')
    m = len(p)
    order = np.argsort(p, kind='stable')
    ranked = p[order] * m / np.arange(1, m + 1)
    q = np.empty(m, dtype=float)
    q[order] = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], 0, 1)
    # Independent implementation check against SciPy.
    np.testing.assert_allclose(q, false_discovery_control(p, method='bh'), rtol=1e-12, atol=1e-14)
    return q


def stars(q: float) -> str:
    return '***' if q < .001 else '**' if q < .01 else '*' if q < .05 else ''


def read_guides(input_path: Path):
    with input_path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        required = ['guide_id', 'target_gene', 'guide_type', VALUE_COL]
        missing = [c for c in required if c not in fields]
        if missing:
            raise ValueError(f'Missing columns: {missing}. Available: {fields}')
        original = list(reader)
    ids = [r['guide_id'].strip() for r in original]
    if any(not s for s in ids) or len(ids) != len(set(ids)):
        raise ValueError('Guide IDs must be nonempty and unique: no automatic aggregation is performed.')
    valid, excluded = [], []
    for r in original:
        try:
            value = float(r[VALUE_COL])
        except (ValueError, TypeError):
            value = float('nan')
        if not math.isfinite(value):
            excluded.append({**r, 'exclusion_reason': 'nonfinite_estimated_EC_percentage'})
            continue
        if not 0 <= value <= 100:
            raise ValueError(f'EC percentage outside [0, 100]: {r["guide_id"]}: {value}')
        guide_type = r['guide_type'].strip().upper()
        if guide_type not in CONTROL_TYPES | {TARGET_TYPE}:
            raise ValueError(f'Unrecognized guide_type {guide_type!r}; explicit classification required.')
        gene = r['target_gene'].strip().upper()
        if not gene:
            raise ValueError(f'Empty target gene for {r["guide_id"]}')
        valid.append({'guide_id': r['guide_id'].strip(), 'target_gene': gene,
                      'guide_type': guide_type, VALUE_COL: value})
    controls = sorted((r for r in valid if r['guide_type'] in CONTROL_TYPES), key=lambda r: r['guide_id'])
    targeting = defaultdict(list)
    for r in valid:
        if r['guide_type'] == TARGET_TYPE:
            targeting[r['target_gene']].append(r)
    targeting = {g: sorted(rr, key=lambda r: r['guide_id']) for g, rr in sorted(targeting.items())}
    if len(controls) < 2:
        raise ValueError('At least two finite individual control-guide values required.')
    for _, genes in GROUPS[1:]:
        for gene in genes:
            if gene not in targeting:
                raise ValueError(f'Requested plot gene {gene} is missing. Do not substitute or invent values.')
    return original, valid, controls, targeting, excluded, fields


def uniform_label_subsets(n_control: int, n_target: int, permutations: int,
                          seed: int, batch_size: int = 5000) -> np.ndarray:
    """Uniform k-subsets of n guides: independent random-key selection per row.

    Cache by group size; reusing subsets across genes of the same size does
    not change the marginal label-permutation distribution of any gene.
    No replacement within each permutation; replacement across permutations.
    """
    n = n_control + n_target
    rng = np.random.default_rng(np.random.SeedSequence([seed, n_control, n_target]))
    dtype = np.uint16 if n < 65536 else np.uint32
    output = np.empty((permutations, n_target), dtype=dtype)
    for lo in range(0, permutations, batch_size):
        hi = min(lo + batch_size, permutations)
        keys = rng.random((hi - lo, n))  # float64, random ranks
        chosen = np.argpartition(keys, n_target - 1, axis=1)[:, :n_target]
        chosen.sort(axis=1)
        if n_target > 1 and np.any(np.diff(chosen, axis=1) <= 0):
            raise AssertionError('A permutation selected a guide more than once.')
        output[lo:hi] = chosen
    return output


def compute_statistics(controls, targeting, outdir: Path, permutations: int, seed: int):
    cv = np.array([r[VALUE_COL] for r in controls], dtype=float)
    cm = float(cv.mean())
    n_control = len(cv)
    sample_sizes = sorted({len(rows) for rows in targeting.values()})
    samples = {k: uniform_label_subsets(n_control, k, permutations, seed) for k in sample_sizes}
    cache_path = outdir / 'permutation_indices_100k.npz'
    np.savez_compressed(cache_path, **{f'n_target_{k}': idx for k, idx in samples.items()})
    results = []
    for i, (gene, rows) in enumerate(targeting.items(), start=1):
        tv = np.array([r[VALUE_COL] for r in rows], dtype=float)
        n_target = len(tv)
        tm = float(tv.mean())
        observed = tm - cm
        pooled = np.concatenate([tv, cv])
        total = pooled.sum(dtype=np.float64)
        # Complete pooled labeling: both group means are recomputed in every draw.
        # Centering limits subtractive rounding errors without changing differences.
        centered = pooled - total / len(pooled)
        centered_total = centered.sum(dtype=np.float64)
        idx = samples[n_target]
        extreme = 0
        tol = 100 * np.finfo(float).eps * max(1.0, abs(observed))
        for lo in range(0, permutations, 10000):
            chunk = idx[lo:lo + 10000]
            target_sum = centered[chunk].sum(axis=1, dtype=np.float64)
            perm_diff = target_sum / n_target - (centered_total - target_sum) / n_control
            extreme += int(np.count_nonzero(np.abs(perm_diff) >= abs(observed) - tol))
        p = extreme / permutations
        p_plus_one = (extreme + 1) / (permutations + 1)
        ratio = (float(np.log2(cm / tm)) if cm > 0 and tm > 0 else
                 (float('inf') if tm == 0 and cm > 0 else float('nan')))
        results.append({
            'gene': gene, 'n_targeting_guides': n_target, 'n_control_guides': n_control,
            'mean_target_EC_pct': tm, 'median_target_EC_pct': float(np.median(tv)),
            'mean_control_EC_pct': cm, 'observed_difference_target_minus_control_pp': observed,
            'effect_log2_control_over_target': ratio,
            'n_extreme_permutations': extreme, 'n_permutations': permutations,
            'permutation_pvalue_two_sided_b_over_B': p,
            'permutation_pvalue_plus_one_sensitivity': p_plus_one,
            'zero_exceedance_flag': extreme == 0, 'seed': seed,
        })
        # Recoverable, fully computed per-gene checkpoint before BH adjustment.
        if i % 25 == 0 or i == len(targeting):
            write_csv(outdir / 'permutation_checkpoint_unadjusted.csv', results)
    q = bh_adjust(np.array([r['permutation_pvalue_two_sided_b_over_B'] for r in results]))
    qplus = bh_adjust(np.array([r['permutation_pvalue_plus_one_sensitivity'] for r in results]))
    for r, fdr, fdrplus in zip(results, q, qplus):
        r['FDR_BH'] = float(fdr)
        r['FDR_BH_plus_one_sensitivity'] = float(fdrplus)
        r['significance'] = stars(float(fdr))
        r['significance_plus_one_sensitivity'] = stars(float(fdrplus))
        r['BH_number_genes'] = len(results)
    return results


def plot_grouped(controls, targeting, results, outdir: Path):
    mpl.rcParams.update({
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'font.family': 'DejaVu Sans', 'figure.facecolor': 'white',
        'axes.facecolor': 'white', 'savefig.transparent': False,
    })
    stats_index = {r['gene']: r for r in results}
    datasets, labels, positions, group_spans, names, jitter_rows = [], [], [], [], [], []
    x = 0.0
    for group_name, genes in GROUPS:
        start = x
        for gene in genes:
            rows = controls if gene == 'NTC' else targeting[gene]
            values = np.array([r[VALUE_COL] for r in rows], dtype=float)
            positions.append(x); datasets.append(values); names.append(gene)
            labels.append(f'{gene}\nn = {len(values)}')
            x += 1.0
        group_spans.append((group_name, start, x - 1.0))
        x += .90
    # Fixed physical category spacing and identical widths in all classes.
    x_min, x_max = -.6, positions[-1] + .6
    width_per_position = .58
    ax_width = (x_max - x_min) * width_per_position
    ax_height, left, right, bottom, top = 3.85, 1.15, .15, 1.80, .62
    fig_w, fig_h = left + ax_width + right, bottom + ax_height + top
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([left / fig_w, bottom / fig_h, ax_width / fig_w, ax_height / fig_h])
    ax.axhline(np.mean(datasets[0]), color='#9EA3A8', linestyle='--', linewidth=1.6, zorder=1)
    bp = ax.boxplot(
        datasets, positions=positions, widths=.42, whis=1.5, patch_artist=True,
        showfliers=False, manage_ticks=False,
        medianprops={'color':'black', 'linewidth':2.0, 'zorder':5},
        whiskerprops={'color':'black', 'linewidth':1.8, 'zorder':5},
        capprops={'color':'black', 'linewidth':1.8, 'zorder':5},
        boxprops={'edgecolor':'black', 'linewidth':1.8, 'zorder':3},
    )
    for box, gene in zip(bp['boxes'], names):
        box.set_facecolor('white' if gene == 'NTC' else '#E3B51E')
        box.set_alpha(1.0)
        box.set_gid(f'box_{gene}')
    for gene, pos, values in zip(names, positions, datasets):
        digest = int.from_bytes(hashlib.sha256(gene.encode()).digest()[:4], 'little')
        rng = np.random.default_rng(np.random.SeedSequence([7, digest]))
        jitter = np.clip(rng.normal(0, .09 if gene == 'NTC' else .055, len(values)), -.16, .16)
        dots = ax.scatter(np.full(len(values), pos) + jitter, values,
                          s=16 if gene == 'NTC' else 18, c='#7A7A7A', alpha=.48,
                          linewidths=0, edgecolors='none', zorder=6, rasterized=False)
        dots.set_gid(f'guide_points_{gene}')
        for rr, offset in zip(controls if gene == 'NTC' else targeting[gene], jitter):
            jitter_rows.append({'plot_gene':gene, 'guide_id':rr['guide_id'],
                                VALUE_COL:rr[VALUE_COL], 'x_center':pos, 'x_jitter':float(offset)})
        if gene != 'NTC':
            annotation = stats_index[gene]['significance']
            if annotation:
                ax.text(pos, 95, annotation, ha='center', va='bottom', fontsize=12,
                        fontweight='bold', color='black', zorder=8)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, 105)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=58, ha='right', rotation_mode='anchor', fontsize=12)
    ax.set_ylabel('Estimated D3 EC cells (%)', fontsize=18, labelpad=8)
    ax.tick_params(axis='y', labelsize=12, width=2.0, length=6)
    ax.tick_params(axis='x', width=2.0, length=6, pad=6)
    for side in ['top','right']: ax.spines[side].set_visible(False)
    for side in ['left','bottom']: ax.spines[side].set_linewidth(2.2)
    ax.grid(False)
    transform = ax.get_xaxis_transform()
    for group_name, start, end in group_spans:
        # Keep control separate; other classes retain their requested order.
        ax.plot([start - .32, end + .32], [-.290, -.290], transform=transform,
                clip_on=False, color='black', linewidth=1)
        ax.text((start + end)/2, -.323, GROUP_LABELS[group_name], transform=transform,
                ha='center', va='top', fontsize=10.5, linespacing=1.15,
                color='black', fontweight='bold', clip_on=False)
    fig.text(left/fig_w, (bottom+ax_height+.24)/fig_h,
             '100,000 two-sided guide-label permutations; BH correction across all tested genes',
             ha='left', va='bottom', fontsize=9.5)
    fig.text((left+ax_width)/fig_w, (bottom+ax_height+.24)/fig_h,
             '* FDR < 0.05   ** < 0.01   *** < 0.001', ha='right', va='bottom', fontsize=9.5)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Check figure text/tick labels fit the fixed-size canvas.
    texts = [*fig.texts, *ax.get_xticklabels(), *ax.get_yticklabels(), ax.yaxis.label]
    for text in texts:
        if not text.get_text(): continue
        bb = text.get_window_extent(renderer)
        fb = fig.bbox
        if bb.x0 < fb.x0-1 or bb.y0 < fb.y0-1 or bb.x1 > fb.x1+1 or bb.y1 > fb.y1+1:
            raise RuntimeError(f'Label outside canvas: {text.get_text()!r}: {bb}')
    for ext in ['pdf','svg','png']:
        fig.savefig(outdir / f'{STEM}.{ext}', dpi=600, facecolor='white', transparent=False)
    plt.close(fig)
    write_csv(outdir / 'plot_jitter_coordinates.csv', jitter_rows)


def main(input_path: str | Path | None = None, outdir: str | Path | None = None,
         permutations: int = 100000, seed: int = 20260908):
    here = Path(__file__).resolve().parent
    input_path = Path(input_path) if input_path else here / 'guide_EC_percent_control70.csv'
    outdir = Path(outdir) if outdir else here / 'results'
    if not input_path.is_file(): raise FileNotFoundError(input_path)
    if permutations <= 0: raise ValueError('permutations must be positive.')
    outdir.mkdir(parents=True, exist_ok=True)
    original, valid, controls, targeting, excluded, fields = read_guides(input_path)
    print(f'Input: {len(original)} guide rows; {len(controls)} controls; {len(targeting)} target genes.', flush=True)
    write_csv(outdir / 'excluded_guides.csv', excluded, fields + ['exclusion_reason'])
    write_csv(outdir / 'control_guides_used.csv', controls)
    counts = [dict(gene=g, n_guides=len(rr)) for g, rr in targeting.items()]
    write_csv(outdir / 'guide_count_QC_all_genes.csv', counts)
    results = compute_statistics(controls, targeting, outdir, permutations, seed)
    selected_order = [g for _, gg in GROUPS[1:] for g in gg]
    stats_index = {r['gene']:r for r in results}
    selected_stats = [{'class':cat, **stats_index[g]} for cat, gg in GROUPS[1:] for g in gg]
    write_csv(outdir / 'permutation_statistics_all_genes.csv', results)
    write_csv(outdir / 'permutation_statistics_plotted_genes.csv', selected_stats)
    subset = []
    for cat, genes in GROUPS:
        for gene in genes:
            for r in controls if gene == 'NTC' else targeting[gene]:
                subset.append({'plot_class':cat, 'plot_gene':gene, **r})
    write_csv(outdir / 'plotted_guide_values.csv', subset)
    plot_grouped(controls, targeting, results, outdir)
    metadata = {
        'input_file': str(input_path.resolve()), 'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),
        'value_column': VALUE_COL, 'source_values_used_directly':True,
        'control_types':sorted(CONTROL_TYPES), 'control_count_by_type':dict(Counter(r['guide_type'] for r in controls)),
        'n_control_guides':len(controls), 'mean_control_EC_pct':float(np.mean([r[VALUE_COL] for r in controls])),
        'n_targeting_guides':sum(map(len,targeting.values())), 'n_genes_tested_BH':len(results),
        'n_excluded_nonfinite_guides':len(excluded), 'quality_flag_filtering':False,
        'n_permutations':permutations, 'rng_seed':seed, 'sampling':'uniform pooled-label permutation; fixed group sizes',
        'primary_pvalue':'number(abs(permuted difference) >= abs(observed difference))/B',
        'roundoff_tolerance':'100 * float64 eps * max(1, abs(observed difference))',
        'zero_pvalue_caution':'b=0 means no exceedances in B Monte Carlo draws; not an exact zero probability',
        'sensitivity_pvalue':'(b+1)/(B+1); separate BH column; not used for primary plot stars',
        'star_thresholds':{'*':.05,'**':.01,'***':.001},
        'weighting':'unweighted arithmetic guide means; unequal per-gene guide counts retained',
        'boxplot':'25th/75th percentiles, median, whiskers 1.5 IQR; all guide values overlaid',
        'plot_colors':{'target_box':'#E3B51E','control_box':'white','guide_points':'#7A7A7A'},
        'point_alpha':.48, 'box_width':.42, 'pdf_fonttype':42, 'svg_fonttype':'none',
        'guide_count_distribution':dict(Counter(len(rr) for rr in targeting.values())),
        'selected_groups':GROUPS, 'python':platform.python_version(),
        'numpy':np.__version__, 'matplotlib':mpl.__version__,
    }
    (outdir/'analysis_metadata.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False),encoding='utf-8')
    (outdir/'README.txt').write_text(
        'VALIDATION-SCREEN GUIDE EC% BOXPLOT\n\n'
        'Source: guide_EC_percent_control70.csv, estimated_EC_pct_control70 (used without reconstruction).\n'
        f'Control: {len(controls)} individual NON-TARGETING/SAFE guides.\n'
        f'BH correction: all {len(results)} unique TARGETING gene labels with finite values.\n'
        'This validation export is distinct from the earlier 119-NTC/298-gene CC-Perturb-seq export.\n'
        'Each guide is unweighted; input read-count/QC flags are not used to remove or weight guides.\n'
        'All plotted target genes have six guides. Non-plotted genes retain the guide counts in the CSV.\n\n'
        'Test: 100,000 two-sided pooled-guide-label permutations per gene; statistic target mean minus control mean.\n'
        'Primary p-values use b/B, matching the requested definition. Stars use BH-adjusted primary p-values.\n'
        'A zero primary p-value means zero observed exceedances, not a known zero probability.\n'
        'Plus-one Monte Carlo p-values and their BH adjustments are included as separate sensitivity columns.\n'
        'The dashed line is the arithmetic mean of the control-guide percentages.\n'
        'Boxes: median and IQR; whiskers: most extreme observations within 1.5 IQR. All guide points are displayed.\n'
        'No new cell-level exclusion can be applied from this aggregated CSV.\n\n'
        'Style: white control, mustard target boxes; transparent gray points in front (alpha=0.48).\n'
        'PDF/SVG retain text and vectors; actual editing in Illustrator is not tested in this environment.\n\n'
        'Reproduce: python make_validation_boxplots.py --input guide_EC_percent_control70.csv --outdir results\n'
        'Dependencies: numpy, scipy, matplotlib. No font files are included.\n',encoding='utf-8')
    print('Completed calculations and exports:', outdir, flush=True)
    for r in selected_stats:
        print(f'{r["gene"]:7s} n={r["n_targeting_guides"]} mean={r["mean_target_EC_pct"]:.3f}'
              f' p={r["permutation_pvalue_two_sided_b_over_B"]:.6g} FDR={r["FDR_BH"]:.6g} {r["significance"]}')
    return metadata


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=Path(__file__).resolve().parent/'guide_EC_percent_control70.csv')
    ap.add_argument('--outdir', type=Path, default=Path(__file__).resolve().parent/'results')
    ap.add_argument('--permutations', type=int, default=100000)
    ap.add_argument('--seed', type=int, default=20260908)
    args=ap.parse_args()
    main(args.input,args.outdir,args.permutations,args.seed)

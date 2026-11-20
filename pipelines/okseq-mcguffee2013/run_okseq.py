#!/usr/bin/env python3
"""Re-analysis of McGuffee et al. 2013 Okazaki-fragment sequencing (GSE40696, wild type).

Produces fork-directionality profiles, initiation-site calls and the OriDB resolution
measurement from mapped reads; validates against the authors' published .sgr signal.

    python run_okseq.py --bam-dir okseq/bam --out results

Inputs
  --bam-dir   directory holding SRR566701.r1.bam / SRR566702.r1.bam (read-1 records of
              properly paired, MAPQ>=20 alignments produced by map.sh)
  --sgr-npz   per-base published Watson/Crick coverage lifted to sacCer3 (validation target)
  --origins   OriDB v3 origins.tsv (Confirmed tier; used only for polarity confirmation,
              the resolution measurement and the circularity check - never for tuning)
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import okseq_lib as K

SAMPLES = {'SRR566701': 'WT_A', 'SRR566702': 'WT_B'}
RDNA = {'chrXII': [(450000, 470000)]}          # rDNA repeat array, sacCer3
PRIMARY = dict(step=250, sigma_kb=2.0, min_sep_bp=5000, oem_quad=10000, oem_cut=0.05)


def profiles(cov, step, sigma_kb, oem_quad=10000):
    win = K.window_sums(cov, step)
    fd = K.directionality(win)                       # (Crick - Watson)/(Crick + Watson)
    oe = K.oem(win, quad_bp=oem_quad, step=step)     # McGuffee OEM = (FD_right - FD_left)/2
    fds = K.smooth_fd(fd, sigma_kb * 1000.0 / step)
    dv = K.derivative(fds, step)
    return win, fd, fds, oe, dv


def call_set(cov, step, sigma_kb, min_sep_bp, oem_quad, oem_cut, mask=RDNA):
    win, fd, fds, oe, dv = profiles(cov, step, sigma_kb, oem_quad)
    c = K.call_sites(fds, dv, step=step, min_height=1e-6, min_sep_bp=min_sep_bp, mask=mask)
    c = K.annotate_oem(c, oe, step)
    c = c[c.oem > oem_cut].reset_index(drop=True)
    return c, (win, fd, fds, oe, dv)


def merge_counts(covs):
    """Pool replicates after equalising each sample's total depth (so neither dominates)."""
    tot = [sum(float(w.sum() + c.sum()) for w, c in cv.values()) for cv in covs]
    ref = float(np.mean(tot))
    out = {}
    for c in K.CHROMS:
        W = sum(cv[c][0] * (ref / t) for cv, t in zip(covs, tot))
        C = sum(cv[c][1] * (ref / t) for cv, t in zip(covs, tot))
        out[c] = (W.astype(np.float32), C.astype(np.float32))
    return out


def correlate(a, b, chroms=K.CHROMS):
    x = np.concatenate([a[c] for c in chroms]); y = np.concatenate([b[c] for c in chroms])
    m = np.isfinite(x) & np.isfinite(y)
    from scipy.stats import spearmanr
    return float(np.corrcoef(x[m], y[m])[0, 1]), float(spearmanr(x[m], y[m]).statistic), int(m.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam-dir', default='okseq/bam')
    ap.add_argument('--sgr-npz', default='okseq/pub_sgr_sacCer3.npz')
    ap.add_argument('--origins', default='/Users/nieduszc/Documents/oridb-v3/public/data/origins.tsv')
    ap.add_argument('--sizes', default='okseq/ref/sacCer3.chrom.sizes')
    ap.add_argument('--out', default='results')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sizes = K.chrom_sizes(a.sizes)
    ori = pd.read_csv(a.origins, sep='\t')
    conf = ori[(ori.species == 'Saccharomyces cerevisiae') & (ori.status == 'Confirmed')].copy()
    summary = {}

    # ---- coverage from BAM -------------------------------------------------
    covs, stats = {}, {}
    for srr, name in SAMPLES.items():
        bam = os.path.join(a.bam_dir, srr + '.r1.bam')
        cov, nfrag, flen = K.fragment_coverage(bam, sizes)
        covs[name] = cov
        stats[name] = dict(fragments=int(sum(nfrag.values())),
                           frag_len_median=float(np.median(flen)),
                           frag_len_mean=float(flen.mean()))
        print(name, stats[name], flush=True)
    merged = merge_counts([covs['WT_A'], covs['WT_B']])
    summary['fragment_stats'] = stats

    # ---- replicate agreement ----------------------------------------------
    _, fdA, _, oeA, _ = profiles(covs['WT_A'], 1000, 2.0)
    _, fdB, _, oeB, _ = profiles(covs['WT_B'], 1000, 2.0)
    summary['replicate_agreement'] = dict(zip(['pearson', 'spearman', 'n'], correlate(oeA, oeB)))

    # ---- validation against the authors' published signal -------------------
    pub = np.load(a.sgr_npz)
    val = {}
    for s, name in [('A', 'WT_A'), ('B', 'WT_B')]:
        pcov = {c: (pub[f'{s}_Watson_{c}'], pub[f'{s}_Crick_{c}']) for c in K.CHROMS}
        _, pfd, _, poe, _ = profiles(pcov, 1000, 2.0)
        _, mfd, _, moe, _ = profiles(covs[name], 1000, 2.0)
        val[name] = dict(OEM=dict(zip(['pearson', 'spearman', 'n'], correlate(poe, moe))),
                         FD=dict(zip(['pearson', 'spearman', 'n'], correlate(pfd, mfd))))
        # strand-level check (pins the Watson/Crick convention without using origins)
        wa = {c: K.window_sums({c: covs[name][c]}, 1000)[c][0] for c in K.CHROMS}
        pw = {c: K.window_sums({c: pcov[c]}, 1000)[c][0] for c in K.CHROMS}
        pc = {c: K.window_sums({c: pcov[c]}, 1000)[c][1] for c in K.CHROMS}
        val[name]['mine_plus_vs_pub_Watson'] = correlate(wa, pw)[0]
        val[name]['mine_plus_vs_pub_Crick'] = correlate(wa, pc)[0]
    summary['validation_vs_published'] = val
    json.dump(summary, open(os.path.join(a.out, 'summary_part1.json'), 'w'), indent=1)
    print(json.dumps(summary['validation_vs_published'], indent=1), flush=True)
    return covs, merged, conf, sizes, summary


if __name__ == '__main__':
    main()


# ============================ stage 2: scans, calls, tracks ==================

def window_scan(covs, merged, conf, steps=(100, 250, 500, 1000, 2000, 5000),
                sigma_kb=PRIMARY['sigma_kb']):
    """Sensitivity of the result to window size, at the primary smoothing."""
    rows = []
    for st in steps:
        cs = {n: call_set(covs[n], st, sigma_kb, PRIMARY['min_sep_bp'],
                          PRIMARY['oem_quad'], PRIMARY['oem_cut'])[0] for n in covs}
        cm = call_set(merged, st, sigma_kb, PRIMARY['min_sep_bp'],
                      PRIMARY['oem_quad'], PRIMARY['oem_cut'])[0]
        ab, ba = K.concordance(cs['WT_A'], cs['WT_B'])
        d = K.nearest_distance(cm, conf)
        r = K.resolution_summary(d, len(conf))
        rows.append(dict(window_bp=st, n_calls_merged=len(cm), n_WT_A=len(cs['WT_A']),
                         n_WT_B=len(cs['WT_B']), replicate_concordance=round((ab + ba) / 2, 3),
                         median_width_bp=int(cm.width.median()),
                         err_median=round(r['median']), err_p90=round(r['p90']),
                         pct_detected=round(r['pct_detected'], 1)))
    return pd.DataFrame(rows)


def smoothing_scan(covs, merged, conf, sigmas=(0.25, 0.5, 1.0, 2.0, 4.0),
                   step=PRIMARY['step']):
    rows = []
    for sg in sigmas:
        cs = {n: call_set(covs[n], step, sg, PRIMARY['min_sep_bp'],
                          PRIMARY['oem_quad'], PRIMARY['oem_cut'])[0] for n in covs}
        cm = call_set(merged, step, sg, PRIMARY['min_sep_bp'],
                      PRIMARY['oem_quad'], PRIMARY['oem_cut'])[0]
        ab, ba = K.concordance(cs['WT_A'], cs['WT_B'])
        off = matched_offsets(cs['WT_A'], cs['WT_B'])
        d = K.nearest_distance(cm, conf)
        r = K.resolution_summary(d, len(conf))
        rows.append(dict(sigma_kb=sg, n_calls_merged=len(cm),
                         replicate_concordance=round((ab + ba) / 2, 3),
                         xrep_offset_median=int(np.median(np.abs(off))) if len(off) else -1,
                         xrep_offset_p90=int(np.percentile(np.abs(off), 90)) if len(off) else -1,
                         median_width_bp=int(cm.width.median()),
                         err_median=round(r['median']), err_p90=round(r['p90']),
                         pct_detected=round(r['pct_detected'], 1)))
    return pd.DataFrame(rows)


def matched_offsets(a, b, tol=3000):
    offs = []
    for c, g in a.groupby('chrom'):
        h = b[b.chrom == c].peak.values
        if len(h) == 0: continue
        dif = g.peak.values[:, None] - h[None, :]
        j = np.abs(dif).argmin(1)
        dmin = dif[np.arange(len(g)), j]
        offs.append(dmin[np.abs(dmin) <= tol])
    return np.concatenate(offs) if offs else np.array([])


def polarity_evidence(oe_prof, conf, step, seed=0):
    """Is OEM positive at Confirmed origins (origins) or negative (termination zones)?"""
    rng = np.random.default_rng(seed)
    at, rand = [], []
    for c, g in conf.groupby('chromosome'):
        o = oe_prof[c]
        idx = ((g.start + g.end) // 2).values // step
        idx = idx[(idx >= 0) & (idx < len(o))]
        at.append(o[idx]); rand.append(o[rng.integers(0, len(o), size=5 * len(idx))])
    at = np.concatenate(at); rand = np.concatenate(rand)
    at = at[np.isfinite(at)]; rand = rand[np.isfinite(rand)]
    return dict(mean_oem_at_confirmed=float(at.mean()), mean_oem_random=float(rand.mean()),
                frac_positive_at_confirmed=float((at > 0).mean()),
                frac_positive_random=float((rand > 0).mean()), n=int(len(at)))


def circularity_check(calls, conf, origins_all):
    """Could these calls have been curated onto the OriDB list? (cf. Crabbe2010)"""
    ident = midpt = 0
    for _, o in conf.iterrows():
        g = calls[calls.chrom == o.chromosome]
        if len(g) == 0: continue
        if ((g.start.values == o.start) & (g.end.values == o.end)).any(): ident += 1
        if (g.peak.values == (o.start + o.end) // 2).any(): midpt += 1
    return dict(n_calls=len(calls), n_confirmed=len(conf),
                pct_calls_interval_identical=round(100.0 * ident / max(1, len(calls)), 3),
                pct_calls_midpoint_identical=round(100.0 * midpt / max(1, len(calls)), 3))


def write_bed(calls, path):
    b = calls.copy()
    b['name'] = ['McGuffee2013_init_%04d' % (i + 1) for i in range(len(b))]
    b['score'] = np.clip((b.oem * 1000).round(), 0, 1000).astype(int)
    cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'peak_start', 'peak_end',
            'itemRgb', 'oem', 'delta_fd', 'deriv']
    b['strand'] = '.'; b['peak_start'] = b.peak; b['peak_end'] = b.peak + 1; b['itemRgb'] = '0,0,0'
    with open(path, 'w') as fh:
        fh.write('track name=McGuffee2013_Okseq_initiation description="Ok-seq initiation '
                 'zones, McGuffee 2013 GSE40696 wild type, re-analysed on sacCer3" useScore=1\n')
        fh.write('#' + '\t'.join(cols) + '\n')
        b[cols].to_csv(fh, sep='\t', header=False, index=False)
    return b

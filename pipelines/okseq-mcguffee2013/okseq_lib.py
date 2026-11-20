
"""Ok-seq analysis: fragment coverage per strand, fork directionality, OEM, origin calls."""
import numpy as np, pandas as pd, subprocess, io, os, shutil
SAMTOOLS = shutil.which('samtools') or '/Users/nieduszc/.claude-science/conda/envs/okseq/bin/samtools'

ROMAN = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
CHROMS = ['chr'+r for r in ROMAN]

def chrom_sizes(path):
    return {k: int(v) for k, v in (l.split() for l in open(path).read().strip().split('\n'))}

def fragment_coverage(bam, sizes, chroms=CHROMS, mapq=20):
    """Per-base fragment coverage split by the strand of read 1.
    Returns dict chrom -> (cov_plus, cov_minus) float32 arrays, plus a fragment count/length table."""
    cov, frag_n, frag_len = {}, {}, []
    for c in chroms:
        n = sizes[c]
        p = subprocess.Popen([SAMTOOLS, 'view', '-@', '2', '-q', str(mapq), bam, c],
                             stdout=subprocess.PIPE)
        d_plus = np.zeros(n + 2, dtype=np.int32); d_minus = np.zeros(n + 2, dtype=np.int32)
        nfrag = 0; lens = []
        for chunk in pd.read_csv(p.stdout, sep='\t', header=None, usecols=[1, 3, 7, 8],
                                 names=['flag', 'pos', 'pnext', 'tlen'],
                                 dtype={'flag': np.int32, 'pos': np.int64, 'pnext': np.int64,
                                        'tlen': np.int64}, chunksize=4_000_000, quoting=3):
            flag = chunk['flag'].values; pos = chunk['pos'].values
            pnext = chunk['pnext'].values; tlen = chunk['tlen'].values
            left = np.where(tlen > 0, pos, pnext) - 1          # 0-based fragment start
            flen = np.abs(tlen)
            ok = (flen > 0) & (left >= 0) & (left + flen <= n)
            left = left[ok]; flen = flen[ok]; rev = (flag[ok] & 16) > 0
            nfrag += ok.sum(); lens.append(flen[:200000])
            for sel, darr in ((~rev, d_plus), (rev, d_minus)):
                s = left[sel]; e = s + flen[sel]
                np.add.at(darr, s, 1); np.add.at(darr, e, -1)
        p.stdout.close(); p.wait()
        cov[c] = (np.cumsum(d_plus)[:n].astype(np.float32), np.cumsum(d_minus)[:n].astype(np.float32))
        frag_n[c] = int(nfrag)
        if lens: frag_len.append(np.concatenate(lens))
    return cov, frag_n, (np.concatenate(frag_len) if frag_len else np.array([]))

def window_sums(cov, step=1000):
    """Sum per-base coverage into fixed windows. Returns dict chrom -> (W, C) arrays per window."""
    out = {}
    for c, (a, b) in cov.items():
        n = len(a) // step * step
        out[c] = (a[:n].reshape(-1, step).sum(1), b[:n].reshape(-1, step).sum(1))
        del n
    return out

def directionality(win, watson_is_plus=True, pseudo=0.0):
    """Fork directionality FD = (C - W)/(C + W) per window.
    W = Watson-strand (leftward-fork) coverage, C = Crick-strand (rightward-fork) coverage.
    FD = -1 where all forks move leftward, +1 where all move rightward."""
    fd = {}
    for c, (p, m) in win.items():
        W, C = (p, m) if watson_is_plus else (m, p)
        tot = W + C
        with np.errstate(invalid='ignore', divide='ignore'):
            f = (C - W) / np.where(tot > 0, tot, np.nan)
        fd[c] = f.astype(np.float32)
    return fd

def oem(win, quad_bp=10000, step=1000, watson_is_plus=True):
    """Origin Efficiency Metric of McGuffee 2013: OEM = WLn - WRn,
    WLn = WL/(WL+CL) over the quadrant left of x, WRn = WR/(WR+CR) over the quadrant right.
    Positive at origins (leftward forks to the left, rightward forks to the right)."""
    q = quad_bp // step
    out = {}
    for c, (p, m) in win.items():
        W, C = (p, m) if watson_is_plus else (m, p)
        cW = np.concatenate([[0.0], np.cumsum(W, dtype=np.float64)])
        cC = np.concatenate([[0.0], np.cumsum(C, dtype=np.float64)])
        n = len(W)
        res = np.full(n, np.nan)
        i = np.arange(q, n - q + 1)                      # boundary index between quadrants
        WL = cW[i] - cW[i - q]; CL = cC[i] - cC[i - q]
        WR = cW[i + q] - cW[i]; CR = cC[i + q] - cC[i]
        with np.errstate(invalid='ignore', divide='ignore'):
            WLn = WL / np.where(WL + CL > 0, WL + CL, np.nan)
            WRn = WR / np.where(WR + CR > 0, WR + CR, np.nan)
        res[i - 1] = WLn - WRn                            # value assigned to window left of boundary
        out[c] = res.astype(np.float32)
    return out


# ---------------------------------------------------------------- calling ----
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths

def smooth_fd(fd, sigma_win):
    """Gaussian smoothing of the FD profile; sigma_win in units of windows. NaNs interpolated."""
    out = {}
    for c, f in fd.items():
        x = f.astype(np.float64).copy()
        bad = ~np.isfinite(x)
        if bad.all():
            out[c] = x; continue
        idx = np.arange(len(x))
        x[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
        out[c] = gaussian_filter1d(x, sigma_win, mode='nearest') if sigma_win > 0 else x
    return out

def derivative(fds, step=1000):
    """d(FD)/dx, per kb."""
    return {c: np.gradient(f) * (1000.0 / step) for c, f in fds.items()}

def call_sites(fds, deriv, step=1000, min_height=0.02, min_sep_bp=5000, mask=None):
    """Peaks in d(FD)/dx = initiation sites. Width = FWHM of the derivative peak.
    Strength = (FD at right FWHM edge - FD at left edge)/2, the local efficiency estimate."""
    rows = []
    for c in sorted(deriv, key=lambda x: CHROMS.index(x)):
        d = deriv[c].copy(); f = fds[c]
        if mask and c in mask:
            for a, b in mask[c]:
                d[max(0, a // step):min(len(d), b // step + 1)] = np.nan
        dd = np.where(np.isfinite(d), d, -np.inf)
        pk, props = find_peaks(dd, height=min_height, distance=max(1, min_sep_bp // step))
        if len(pk) == 0: continue
        floor = np.nanmin(d[np.isfinite(d)]) if np.isfinite(d).any() else 0.0
        dw = np.where(np.isfinite(d), d, floor)          # finite copy for width measurement
        w, wh, lf, rt = peak_widths(dw, pk, rel_height=0.5)
        lf = np.where(np.isfinite(lf), lf, pk - 1); rt = np.where(np.isfinite(rt), rt, pk + 1)
        for i, p in enumerate(pk):
            l = int(np.floor(lf[i])); r = int(np.ceil(rt[i]))
            l = max(0, l); r = min(len(f) - 1, r)
            rows.append(dict(chrom=c, start=l * step, end=(r + 1) * step,
                             peak=int(p * step + step // 2), width=(r - l + 1) * step,
                             deriv=float(dd[p]),
                             delta_fd=float(f[r] - f[l]),
                             strength=float((f[r] - f[l]) / 2.0)))
    return pd.DataFrame(rows)

def nearest_distance(calls, origins, window=20000, mode='midpoint', pos_col='peak'):
    """Distance from each Confirmed origin to that study's nearest call, within `window`.
    mode='midpoint': |origin midpoint - call position|  (the definition used in
    oridb-v3/docs/study_resolution.csv, where width and error are separate axes)
    mode='gap':      interval-to-interval gap, 0 when the call interval covers the origin."""
    out = []
    for _, o in origins.iterrows():
        g = calls[calls.chrom == o.chromosome]
        if len(g) == 0:
            out.append(np.nan); continue
        omid = (o.start + o.end) // 2
        if mode == 'midpoint':
            d = np.abs(g[pos_col].values - omid)
        else:
            d = np.maximum(0, np.maximum(g.start.values - o.end, o.start - g.end.values))
        m = float(d.min())
        out.append(m if m <= window else np.nan)
    return np.array(out, dtype=float)

def resolution_summary(d, n_total):
    ok = d[np.isfinite(d)]
    return dict(n_detected=int(len(ok)), pct_detected=100.0 * len(ok) / n_total,
                median=float(np.median(ok)), p75=float(np.percentile(ok, 75)),
                p90=float(np.percentile(ok, 90)), p95=float(np.percentile(ok, 95)))

def write_bedgraph(prof, path, step=1000, sizes=None):
    with open(path, 'w') as fh:
        for c in sorted(prof, key=lambda x: CHROMS.index(x)):
            v = prof[c]
            for i in np.flatnonzero(np.isfinite(v)):
                s = i * step; e = min(s + step, sizes[c]) if sizes else s + step
                if e > s: fh.write(f"{c}\t{s}\t{e}\t{v[i]:.5f}\n")



def write_bigwig(prof, path, sizes, step=250, chroms=CHROMS):
    """Write a per-window profile dict (chrom -> array) to bigWig.

    One interval per window; windows whose value is not finite are omitted rather than
    written as zero, so a gap in the data stays a gap in the browser. Written with
    pyBigWig because bedGraphToBigWig in bioconda is broken on this platform (missing
    dylib); the format is the same either way.
    """
    import pyBigWig
    bw = pyBigWig.open(path, 'w')
    bw.addHeader([(c, sizes[c]) for c in chroms])
    for c in chroms:
        v = np.asarray(prof[c], dtype=np.float64)
        idx = np.flatnonzero(np.isfinite(v))
        if len(idx) == 0:
            continue
        starts = idx * step
        ends = np.minimum(starts + step, sizes[c])
        keep = ends > starts
        starts, ends, vals = starts[keep], ends[keep], v[idx][keep]
        bw.addEntries([c] * len(starts), starts.tolist(),
                      ends=ends.tolist(), values=vals.tolist())
    bw.close()
    return path

def annotate_oem(calls, oem_prof, step=1000):
    """Attach the McGuffee OEM (10 kb quadrants) at each call position as the efficiency estimate."""
    v = []
    for _, r in calls.iterrows():
        o = oem_prof[r.chrom]; i = int(r.peak // step)
        v.append(float(o[i]) if 0 <= i < len(o) and np.isfinite(o[i]) else np.nan)
    out = calls.copy(); out['oem'] = v
    return out

def concordance(a, b, tol=2000):
    """Fraction of calls in each set with a partner in the other within tol bp."""
    def frac(x, y):
        n = 0
        for c, g in x.groupby('chrom'):
            h = y[y.chrom == c].peak.values
            if len(h) == 0: continue
            n += int((np.abs(g.peak.values[:, None] - h[None, :]).min(1) <= tol).sum())
        return n / max(1, len(x))
    return frac(a, b), frac(b, a)

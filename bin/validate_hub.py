#!/usr/bin/env python3
"""Check a track hub is internally consistent before it is served.

A broken hub fails in the browser with a message that names a line number in a file the
reader has never seen, so it is worth catching these here: a bigDataUrl pointing at a file
that is not there, a track whose sequence names do not match the assembly, a subtrack whose
parent does not exist, or a bigWig whose values fall outside the declared viewLimits.
"""
import argparse, os, re, sys

BIGWIG_MAGIC = 0x888FFC26
TOTAL_SUMMARY_OFFSET_POS = 44        # bbiFile header: where totalSummaryOffset lives


def bigwig_summary(path):
    """Read a bigWig's total-summary block directly.

    pyBigWig 0.3.30's header() returns minVal and maxVal cast to int, so a profile
    spanning 0.933 to 2.276 is reported as 0 to 2. The file itself is correct - only
    that accessor is lossy - so anything that checks a track's value range has to read
    the summary block rather than trust header().
    """
    import struct
    with open(path, "rb") as fh:
        head = fh.read(64)
        if struct.unpack("<I", head[:4])[0] == BIGWIG_MAGIC:
            endian = "<"
        elif struct.unpack(">I", head[:4])[0] == BIGWIG_MAGIC:
            endian = ">"
        else:
            raise ValueError(f"{path} is not a bigWig file")
        offset = struct.unpack(endian + "Q",
                               head[TOTAL_SUMMARY_OFFSET_POS:TOTAL_SUMMARY_OFFSET_POS + 8])[0]
        if offset == 0:
            raise ValueError(f"{path} has no total-summary block")
        fh.seek(offset)
        valid, vmin, vmax, total, total_sq = struct.unpack(endian + "Qdddd", fh.read(40))
    return {"validCount": valid, "minVal": vmin, "maxVal": vmax,
            "sumData": total, "sumSquared": total_sq,
            "mean": total / valid if valid else float("nan")}


def parse_trackdb(path, root):
    """Return stanzas as dicts, following `include` lines."""
    stanzas, cur = [], {}
    for raw in open(path):
        line = raw.strip()
        if line.startswith("include "):
            inc = os.path.join(os.path.dirname(path), line.split(None, 1)[1])
            stanzas += parse_trackdb(inc, root)
            continue
        if not line or line.startswith("#"):
            if cur:
                stanzas.append(cur); cur = {}
            continue
        key, _, value = line.partition(" ")
        if key == "track" and cur:
            stanzas.append(cur); cur = {}
        cur[key] = value.strip()
        cur.setdefault("_dir", os.path.dirname(path))
    if cur:
        stanzas.append(cur)
    return stanzas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hub-root", required=True)
    ap.add_argument("--chrom-sizes", help="assembly chrom.sizes to check track sequences against")
    a = ap.parse_args()
    errors, warnings = [], []

    for required in ("hub.txt", "genomes.txt"):
        if not os.path.exists(os.path.join(a.hub_root, required)):
            errors.append(f"missing {required}")
    if errors:
        print("\n".join(f"ERROR  {e}" for e in errors)); sys.exit(1)

    genomes = []
    for line in open(os.path.join(a.hub_root, "genomes.txt")):
        if line.startswith("genome "):
            genomes.append(line.split()[1].strip())

    sizes = {}
    if a.chrom_sizes:
        for line in open(a.chrom_sizes):
            c, n = line.split()[:2]
            sizes[c] = int(n)

    total = 0
    for genome in genomes:
        tdb = os.path.join(a.hub_root, genome, "trackDb.txt")
        if not os.path.exists(tdb):
            errors.append(f"{genome}: no trackDb.txt"); continue
        stanzas = parse_trackdb(tdb, a.hub_root)
        names = {s["track"] for s in stanzas if "track" in s}

        for s in stanzas:
            name = s.get("track", "<unnamed>")
            if "parent" in s:
                parent = s["parent"].split()[0]
                if parent not in names:
                    errors.append(f"{name}: parent {parent!r} is not defined")
            if "bigDataUrl" in s:
                total += 1
                # Resolve the way UCSC does: relative to the trackDb declared in
                # genomes.txt, not relative to the file an `include` pulled the stanza
                # from. Getting this wrong makes a hub validate locally and then fail in
                # the browser with "couldn't open <genome>/<file>".
                p = os.path.join(os.path.dirname(tdb), s["bigDataUrl"])
                if not os.path.exists(p):
                    errors.append(f"{name}: bigDataUrl {s['bigDataUrl']} resolves to "
                                  f"{os.path.relpath(p, a.hub_root)}, which does not exist "
                                  f"(paths are relative to {genome}/trackDb.txt)")
                    continue
                try:
                    import pyBigWig
                except ImportError:
                    warnings.append("pyBigWig not installed: bigWig contents not checked")
                    continue
                bw = pyBigWig.open(p)
                chroms = bw.chroms()
                bw.close()
                hdr = bigwig_summary(p)   # not bw.header(): see bigwig_summary
                if sizes:
                    bad = {c: n for c, n in chroms.items()
                           if c not in sizes or sizes[c] != n}
                    if bad:
                        errors.append(f"{name}: sequences absent from the assembly or of the "
                                      f"wrong length: {sorted(bad)[:4]}")
                vl = s.get("viewLimits") or _inherited(stanzas, s, "viewLimits")
                negate = (s.get("negateValues") or
                          _inherited(stanzas, s, "negateValues")) == "on"
                if negate:
                    # The axis is the negated data, so compare against negated values -
                    # and negating swaps which end is the minimum.
                    hdr = dict(hdr, minVal=-hdr["maxVal"], maxVal=-hdr["minVal"],
                               mean=-hdr["mean"])
                if vl and ":" in vl:
                    lo, hi = (float(x) for x in vl.split(":"))
                    # Some clipping is intentional: a smoothed copy-number profile
                    # overshoots the 1-2 scale slightly, and 1-2 is the scale that means
                    # something biologically. What is not intentional is limits that cut
                    # off a real part of the distribution, which shows up as the track's
                    # MEAN sitting outside them, or as a large fraction of the range lost.
                    span = hi - lo
                    over = max(0.0, hdr["maxVal"] - hi) + max(0.0, lo - hdr["minVal"])
                    if not (lo <= hdr["mean"] <= hi):
                        errors.append(f"{name}: mean {hdr['mean']:.2f} lies outside "
                                      f"viewLimits {vl} - the limits are wrong for this track")
                    elif over > 0.25 * span:
                        warnings.append(f"{name}: viewLimits {vl} clip {over:.2f} of a "
                                        f"{hdr['minVal']:.2f}-{hdr['maxVal']:.2f} range; "
                                        f"widen them unless the clipping is deliberate")

    for w in warnings:
        print(f"WARN   {w}")
    for e in errors:
        print(f"ERROR  {e}")
    if errors:
        sys.exit(1)
    print(f"OK: {total} tracks across {len(genomes)} assembly(ies), no inconsistencies")


def _inherited(stanzas, child, key):
    parent = child.get("parent", "").split()[0] if "parent" in child else None
    for s in stanzas:
        if s.get("track") == parent:
            return s.get(key)
    return None


if __name__ == "__main__":
    main()

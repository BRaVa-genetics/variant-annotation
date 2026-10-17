#!/usr/bin/env python3
"""Positive control for GPU scoring: re-score variants the precomputed release already covers and compare.

pick:    300 rows of your own lookup output (lookup_chr.sh: Illumina's precomputed scores for your variants), half of
         them at a score >= 0.2 so the cut-offs are tested; written with a header like a score.sh input.
compare: our scores against those precomputed ones. Not bit-identical by design (fp16 and torch.compile; the model
         README says so), so the test is a tolerance: PASS if >= 99% agree within 0.02 and, at 0.10 and 0.40 (the
         cut-offs brava_create_annot.py uses), every variant more than 0.02 from the cut-off gets the same call.
         Precomputed scores are rounded to 2 decimals, so calls within 0.02 of a cut-off may flip: reported only.
         On 3,592 chr21 variants (2026-10-12): median difference 0, max 0.01, no clear-cut call differed.

usage: control.py pick SCORES_DIR OUT.tsv [N=300] [SEED=1]
       control.py compare EXPECTED.tsv OURS.scores.tsv.gz
"""
import glob, gzip, os, random, sys


def pick(scores_dir, out, n="300", seed="1"):
    n = int(n); rng = random.Random(int(seed)); hi, lo = [], []
    for p in sorted(glob.glob(os.path.join(scores_dir, "chr*.tsv.gz"))):
        with gzip.open(p, "rt") as fh:
            for l in fh:
                f = l.rstrip("\n").split("\t")
                if len(f) == 6 and f[1].isdigit():
                    (hi if float(f[5]) >= 0.2 else lo).append(f)
    if len(hi) + len(lo) < n:
        sys.exit("control.py: only %d precomputed rows in %s; run lookup_chr.sh first" % (len(hi) + len(lo), scores_dir))
    k = min(len(hi), n // 2)
    sel = rng.sample(hi, k) + rng.sample(lo, n - k)
    with open(out, "w") as fo:
        fo.write("chrom\tpos\tref\talt\tstrand\tspliceai2_summary_score\n")
        fo.writelines("\t".join(f) + "\n" for f in sel)
    print("control.py: %d control variants (%d with precomputed >= 0.2) -> %s" % (n, k, out))


def compare(expected, ours_path):
    exp = {}
    for l in open(expected):
        f = l.rstrip("\n").split("\t")
        if f[0] != "chrom":
            exp[tuple(f[:5])] = float(f[5])
    ours = {}
    with gzip.open(ours_path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t"); iS = h.index("spliceai2_summary_score")
        for l in fh:
            f = l.rstrip("\n").split("\t"); ours[tuple(f[:5])] = float(f[iS])
    common = [k for k in exp if k in ours]; n = len(common)
    if n == 0:
        print("CONTROL FAIL: no control variant was scored"); return False
    d = sorted(abs(exp[k] - ours[k]) for k in common)
    print("controls %d, scored %d; |ours - precomputed| median %.4f, max %.4f; within 0.02: %.2f%%" % (
        len(exp), n, d[n // 2], d[-1], 100.0 * sum(x <= 0.02 for x in d) / n))
    ok = n == len(exp) and sum(x <= 0.02 for x in d) >= 0.99 * n
    for t in (0.10, 0.40):
        far = [k for k in common if abs(exp[k] - t) > 0.02]
        bad = sum((exp[k] >= t) != (ours[k] >= t) for k in far)
        print("  at %.2f: %d of %d clear-cut calls differ" % (t, bad, len(far)))
        ok = ok and bad == 0
    print("CONTROL %s" % ("PASS" if ok else "FAIL"))
    return ok


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("pick", "compare"):
        sys.exit(__doc__)
    if sys.argv[1] == "pick":
        pick(*sys.argv[2:])
    else:
        sys.exit(0 if compare(*sys.argv[2:4]) else 1)

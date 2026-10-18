"""AlphaMissense-threshold mask: a SAIGE group file holding only missense variants whose AlphaMissense
score reaches the ACMG/AMP PP3 Moderate threshold, 0.906 (Bergquist et al. 2025; on AoU v8 it carries
the same local evidence as BRaVa's REVEL >= 0.773: LR+ 15.4 vs 15.3, genebass-burden
controls/AM_ACMG_CUTOFF.md). A separate file because a group file gives a variant one label per gene,
and every one of these already carries damaging_missense, other_missense or the SpliceAI label there.

Scores come from the FlexRV builder's --provenance output for the same chromosome (alphamissense_files.sh runs
both), so this mask and the flexrv_AM file use one score per (gene, variant). Kept: SOURCE is a real AlphaMissense
row (not 'missing' = gene-mean placeholder, not 'dropped_nonmissense'; with --nonMissense drop the CSQ includes
missense_variant), AM >= threshold. Variant order and IDs follow the group file brava_create_annot.py wrote.

usage: am_mask.py GROUP_FILE.txt.gz PROVENANCE.tsv OUT.txt.gz [THRESHOLD=0.906] [LABEL=am_moderate]
"""
import gzip, sys
from collections import Counter


def main(group, prov, out, thr="0.906", label="am_moderate"):
    thr = float(thr); keep = {}; n = Counter()
    with open(prov) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for l in fh:
            f = l.rstrip("\n").split("\t")
            src, am = f[ix["SOURCE"]], f[ix["AM"]]
            n["missense-bin rows"] += 1
            if src in ("missing", "dropped_nonmissense") or am == "NA":
                n["no AlphaMissense row (" + src + ")"] += 1; continue
            # missense when missense_variant is one of the '&'-joined VEP terms: Tier 1 rows write them all
            # (missense_variant&splice_region_variant), AoU rows only the most severe (fix 2026-10-15)
            if "missense_variant" not in f[ix["CSQ"]].split("&"):
                n["CSQ not missense_variant"] += 1; continue
            if float(am) >= thr:
                keep[(f[ix["GENE"]], f[ix["ID"]])] = 1
    op = gzip.open if out.endswith(".gz") else open
    with gzip.open(group, "rt") as fg, op(out, "wt") as fo:
        var = None
        for l in fg:
            t = l.split()
            if len(t) > 1 and t[1] == "var":
                var = (t[0], t[2:]); continue
            if len(t) > 1 and t[1] == "anno" and var and var[0] == t[0]:
                vs = [v for v in var[1] if (t[0], v) in keep]
                if vs:
                    fo.write("%s var %s\n%s anno %s\n" % (t[0], " ".join(vs), t[0], " ".join([label] * len(vs))))
                    n["regions written"] += 1; n["variants written"] += len(vs)
    if n["variants written"] != len(keep):
        sys.exit("am_mask: %d variants pass the threshold but %d were written -- provenance and group file "
                 "disagree (built from different inputs?)" % (len(keep), n["variants written"]))
    for k, v in sorted(n.items()):
        print("  %10d  %s" % (v, k))
    print("am_mask: AM >= %g -> %d variants in %d regions (%s)" % (thr, n["variants written"], n["regions written"], out))


if __name__ == "__main__":
    main(*sys.argv[1:])

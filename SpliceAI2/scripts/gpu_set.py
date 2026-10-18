#!/usr/bin/env python3
"""The keys SpliceAI2's precomputed release lacks that are worth scoring on a GPU, as input chunks for score.sh.

Kept: missing keys (lookup_chr.sh's CHROM.KIND.missing.tsv) inside a transcript SpliceAI v1 scores
(data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz, this repository) on the key's own strand, AND within
W bp of a GENCODE v39 exon boundary. The release covers all SNVs and ~150M population indels in protein-coding
genes, so what is missing is mostly rare indels; the window keeps exome-capture variants (on chr21 94% of SpliceAI2
>= 0.25 lie within 250 bp). Everything else without a SpliceAI2 score has no splice route.
Rows are sorted by chrom, pos, strand and a chunk never splits a (chrom, pos, strand) group.

usage: gpu_set.py SCORES_DIR GENCODE_V39.gtf.gz V1_ANNOTATION.txt.gz OUT_DIR [W=250] [CHUNK=70000]
"""
import bisect, glob, gzip, os, re, sys
from collections import Counter, defaultdict


def merged_spans(v1_annotation):
    span = defaultdict(list)
    with gzip.open(v1_annotation, "rt") as fh:
        for l in fh:
            if l[0] == "#":
                continue
            f = l.split("\t")
            span[(f[1], f[2])].append((int(f[3]) + 1, int(f[4])))   # TX_START is 0-based in SpliceAI's format
    out = {}
    for k, iv in span.items():
        iv.sort(); m = []
        for s, e in iv:
            if m and s <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[k] = ([s for s, _ in m], [e for _, e in m])
    return out


def main(scores_dir, gtf, v1_annotation, out_dir, W="250", CH="70000"):
    W, CH = int(W), int(CH)
    merged = merged_spans(v1_annotation)

    def inside(c, p, st):
        if (c, st) not in merged:
            return False
        S, E = merged[(c, st)]; i = bisect.bisect_right(S, p) - 1
        return i >= 0 and p <= E[i]
    ex = defaultdict(set)
    with gzip.open(gtf, "rt") as fh:
        for l in fh:
            if l[0] != "#":
                f = l.split("\t", 5)
                if f[2] == "exon":
                    ex[f[0]].update((int(f[3]), int(f[4])))
    ex = {c: sorted(v) for c, v in ex.items()}

    def dist(c, p):
        b = ex.get(c, []); i = bisect.bisect_left(b, p)
        return min((abs(p - b[j]) for j in (i - 1, i) if 0 <= j < len(b)), default=10 ** 9)
    order = ["chr%s" % x for x in list(range(1, 23)) + ["X", "Y"]]
    rows = []; n = Counter()
    for p in sorted(glob.glob(os.path.join(scores_dir, "chr*.missing.tsv"))):
        c, kind = re.match(r"(chr[^.]+)\.(snv|indel)\.missing\.tsv$", os.path.basename(p)).groups()
        for l in open(p):
            pos, ref, alt, st = l.rstrip("\n").split("\t")
            n[kind + " missing"] += 1
            if inside(c, int(pos), st) and dist(c, int(pos)) <= W:
                rows.append((order.index(c), int(pos), st, ref, alt, c)); n[kind + " to score"] += 1
    rows.sort()
    os.makedirs(out_dir, exist_ok=True)
    i = k = 0
    while i < len(rows):
        j = min(i + CH, len(rows))
        while j < len(rows) and rows[j][:3] == rows[j - 1][:3]:
            j += 1
        with open(os.path.join(out_dir, "chunk_%03d.tsv" % k), "w") as fo:
            fo.write("chrom\tpos\tref\talt\tstrand\n")
            for r in rows[i:j]:
                fo.write("%s\t%d\t%s\t%s\t%s\n" % (r[5], r[1], r[3], r[4], r[2]))
        i = j; k += 1
    print(dict(n))
    print("%d keys to score in %d chunk(s) in %s" % (len(rows), k, out_dir))


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])

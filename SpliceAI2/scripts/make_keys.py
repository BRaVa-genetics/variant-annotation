#!/usr/bin/env python3
"""SpliceAI2 lookup keys: every (variant, strand) your annotation needs, per chromosome.

SpliceAI2's precomputed scores (Hugging Face illumina-ai/SpliceAI2-data, precomputed_scores_v2.0) are keyed by
chrom, pos, ref, alt, strand: one score per strand, no gene. An entry (gene, variant) takes the score on its gene's
GENCODE v39 strand, so the keys are the (variant, gene strand) pairs of the entries brava_create_annot.py keeps:
MANE Select or canonical, protein-coding rows of the VEP table (step 1).

Writes OUT_DIR/chrN.snv.keys.tsv and chrN.indel.keys.tsv (pos, ref, alt, strand; sorted by position).
usage: make_keys.py GENCODE_V39.gtf.gz OUT_DIR VEP_TABLE [VEP_TABLE ...]
"""
import gzip, os, re, sys
from collections import Counter, defaultdict


def opn(p):
    return gzip.open(p, "rt") if p.endswith((".gz", ".bgz")) else open(p)


def main(gtf, out_dir, tables):
    strand = {}
    with opn(gtf) as fh:
        for l in fh:
            if l[0] != "#":
                f = l.split("\t", 8)
                if f[2] == "gene":
                    strand[re.search(r'gene_id "([^"]+)"', f[8]).group(1).split(".")[0]] = f[6]
    keys = defaultdict(set); n = Counter()
    for t in tables:
        with opn(t) as fh:
            h = fh.readline().split(); ix = {k: i for i, k in enumerate(h)}
            for l in fh:
                f = l.split()
                if not f:
                    continue
                mane = f[ix["MANE_SELECT"]] not in (".", "NA", "")
                if f[ix["BIOTYPE"]] != "protein_coding" or not (mane or f[ix["CANONICAL"]] == "YES"):
                    continue
                s = strand.get(f[ix["GENE"]].split(".")[0])
                if s is None:
                    n["entries whose gene is not in the GTF"] += 1; continue
                c, p, r, a = re.split(r"[:_]", f[ix["SNP_ID"]])
                c = c if c.startswith("chr") else "chr" + c
                keys[c].add((int(p), r, a, s)); n["entries"] += 1
    os.makedirs(out_dir, exist_ok=True)
    for c, ks in sorted(keys.items()):
        for kind in ("snv", "indel"):
            sel = sorted(k for k in ks if (len(k[1]) == 1 and len(k[2]) == 1) == (kind == "snv"))
            with open(os.path.join(out_dir, "%s.%s.keys.tsv" % (c, kind)), "w") as fo:
                fo.writelines("%d\t%s\t%s\t%s\n" % k for k in sel)
            print("%s %s: %d keys" % (c, kind, len(sel)))
    print(dict(n))


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3:])

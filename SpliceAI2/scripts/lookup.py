"""Join one chromosome of SpliceAI2 precomputed scores (streamed on stdin, already cut to chrom, pos, ref,
alt, strand, spliceai2_summary_score) against your keys (make_keys.py), so the 348 GB release is never
stored: only the matching rows are written, plus OUT.missing.tsv with the keys the release lacks.

stdin columns: chrom pos ref alt strand score (header line allowed).
usage: hf_ranged_cat.sh URL | gzip -dc | cut -f1-5,84 | lookup.py KEYS.tsv OUT.tsv.gz   (lookup_chr.sh does this)
Prints: keys, stream rows read, matched keys, and keys not found (to be scored on a GPU).
"""
import gzip, sys


def main(keys_path, out):
    keys = {}
    with open(keys_path) as fh:
        for l in fh:
            keys[l.rstrip("\n")] = None
    n = 0; hit = 0
    with gzip.open(out + ".tmp", "wt") as fo:
        for l in sys.stdin:
            n += 1
            i = l.index("\t"); j = l.rindex("\t")
            k = l[i + 1:j]
            if k in keys and keys[k] is None:
                keys[k] = l[j + 1:].rstrip("\n"); hit += 1
                fo.write(l)
    import os
    if n == 0:
        # an empty stream is a failed download (seen 2026-10-12: HTTP 401 left an empty output that a
        # resubmission would have skipped as done); keep nothing
        os.remove(out + ".tmp")
        sys.exit("%s: no rows on stdin -- download failed? nothing written" % os.path.basename(out))
    os.replace(out + ".tmp", out)
    miss = [k for k, v in keys.items() if v is None]
    with open(out.replace(".tsv.gz", "") + ".missing.tsv", "w") as fm:
        for k in miss:
            fm.write(k + "\n")
    snv_miss = sum(1 for k in miss if len(k.split("\t")[1]) == 1 and len(k.split("\t")[2]) == 1)
    print("%s: %d keys, %d stream rows, %d matched (%.2f%%), %d missing (%d SNV, %d indel)" % (
        os.path.basename(out), len(keys), n, hit, 100.0 * hit / max(1, len(keys)), len(miss), snv_miss,
        len(miss) - snv_miss), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

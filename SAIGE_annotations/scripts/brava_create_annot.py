#!/usr/bin/env python3
"""BRaVa annotation: VEP 105 + LOFTEE table and SpliceAI2 (or SpliceAI v1) scores -> SAIGE group file + long table.

One label per (gene, variant), on the gene's MANE Select transcript, else its canonical one (protein-coding only).
The first rule that applies wins (October 2026 rules):

  pLoF               LOFTEE high confidence (HC)
  damaging_missense  a missense-type consequence and (REVEL >= 0.773 or CADD >= 28.1); or LOFTEE low confidence (LC)
  spliceai2_gt040    SpliceAI2 >= 0.40 (with --spliceai: spliceai_gt05, SpliceAI v1 DS >= 0.5)
  other_missense     any other missense-type consequence
  synonymous         synonymous with SpliceAI2 < 0.10 (v1: DS < 0.2)
  non_coding         everything else, including synonymous with SpliceAI2 0.10-0.40 (v1: DS 0.2-0.5)

Missense-type: missense_variant, inframe_insertion, inframe_deletion, start_lost, stop_lost,
protein_altering_variant, transcript_amplification. A consequence such as 'missense_variant&splice_region_variant'
belongs to a class when any of its '&'-joined terms does (as BRaVa's script has always done). Non-coding terms are
tested before synonymous ones, so 'synonymous_variant&intron_variant' is non_coding.

What changed from the 2023 rules: SpliceAI no longer routes into damaging_missense (SpliceAI >= 0.2 did); it is its
own category, at a cut-off calibrated to ACMG/AMP PP3 Moderate (the same evidence level as REVEL 0.773 / CADD 28.1);
the labels are the short names; a variant matching no class is non_coding and stays in the group file. A variant
without a SpliceAI score has no splice route. SpliceAI2 and v1 scores are never mixed: one run uses one model.

SpliceAI2 scores each position once per strand, not per gene: an entry takes the score on its gene's strand
(GENCODE v39), and only when the variant lies inside the gene's v39 span (the 'gene-body rule': SpliceAI v1 only
reports genes whose body contains the variant). Otherwise the score is kept in DS_MAX with SPLICEAI_SOURCE
'spliceai2_outside_gene' and gives no splice route.

Outputs:
  OUT                     SAIGE group file: per gene a 'GENE var ...' and a 'GENE anno ...' line (gzip if OUT ends .gz)
  OUT.long.tsv.gz         one row per (gene, variant) with every field its label came from; the --annoTable input of
                          universal-saige's 04_flexrv_groupfile.sh

usage: brava_create_annot.py -v VEP_TABLE (--spliceai2 SCORES.tsv.gz [...] --gtf GENCODE_V39.gtf.gz | -s SPLICEAI.vcf)
                             -w OUT [--cadd_indels CADD.tsv.gz] [--vep_*_col NAME ...]
       brava_create_annot.py --selftest
"""
import argparse, gzip, os, re, sys, tempfile
from collections import Counter, OrderedDict

PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant"]
MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion",
                 "missense_variant", "protein_altering_variant"]
SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]
OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant",
              "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]
INFRAME_CSQS = ["inframe_deletion", "inframe_insertion"]
# LOFTEE gives every one of these an HC or LC call on a protein-coding transcript; a row with one of them and no
# call means LOFTEE failed on that variant (seen with VEP --fork: its SQLite conservation database fails to read)
LOFTEE_CSQS = ["stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"]

REVEL_CUTOFF, CADD_CUTOFF = 0.773, 28.1
# SpliceAI cut-offs (October 2026). SpliceAI2: 0.40 is its ACMG/AMP Moderate cut-off from our calibration (Illumina's
# README pairs 0.25 with v1's 0.5, but 0.25 sits at Supporting); 0.10 is the README's equivalent of v1's 0.2.
# SpliceAI v1: 0.5 (Moderate in the same calibration) and 0.2 (BRaVa's earlier cut-off).
MODELS = {"spliceai2": dict(t_hi=0.40, t_lo=0.10, label="spliceai2_gt040"),
          "spliceai": dict(t_hi=0.50, t_lo=0.20, label="spliceai_gt05")}
NAN = float("nan")
LONG_COLS = ["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED", "CSQ", "TRANSCRIPT", "MANE_SELECT", "CANONICAL",
             "BIOTYPE", "DS_MAX", "SPLICEAI_SOURCE", "ANNOTATION"]
DEFAULT_COLS = dict(vep_snp_id_col="SNP_ID", vep_gene_col="GENE", vep_lof_col="LOF", vep_revel_col="REVEL_SCORE",
                    vep_cadd_phred_col="CADD_PHRED", vep_consequence_col="CSQ", vep_transcript_col="TRANSCRIPT",
                    vep_mane_select_col="MANE_SELECT", vep_canonical_col="CANONICAL", vep_biotype_col="BIOTYPE")


def opn(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith((".gz", ".bgz")) else open(p, mode)


def ge(x, t):
    return x == x and x >= t  # a missing score (NaN) never reaches a cut-off


def classes(csq):
    cs = csq.split("&")
    return (any(c in MISSENSE_CSQS for c in cs), any(c in SYNONYMOUS_CSQS for c in cs),
            any(c in OTHER_CSQS for c in cs), any(c in INFRAME_CSQS for c in cs))


def annotate(csq, lof, revel, cadd, ds, t_hi, t_lo, splice_label):
    """The October 2026 rule. lof: 'HC', 'LC' or ''; revel, cadd, ds: float or NaN."""
    mis, syn, oth, inf = classes(csq)
    if lof == "HC":
        return "pLoF"
    if (mis and (ge(revel, REVEL_CUTOFF) or ge(cadd, CADD_CUTOFF))) or lof == "LC":
        return "damaging_missense"
    if ge(ds, t_hi):
        return splice_label
    if mis or inf:
        return "other_missense"
    if oth:
        return "non_coding"
    if syn:
        return "non_coding" if ge(ds, t_lo) else "synonymous"
    return "non_coding"


def first_value(v):
    """First non-missing number of a list. VEP joins values with '&'; dbNSFP 5.4a stores one REVEL per transcript
    separated by ';' ('.;0.53'), which can also reach the table as ',' or URL-encoded '%3B'."""
    for x in v.replace("%3B", "&").replace(";", "&").replace(",", "&").split("&"):
        if x not in (".", "", "NA"):
            try:
                return float(x)
            except ValueError:
                pass
    return NAN


def norm(vid):
    """chr1:100:A:G and 1:100:A:G (and '_' separators) -> one join key."""
    c, p, r, a = re.split(r"[:_]", vid)
    return "%s:%s:%s:%s" % (c[3:] if c.startswith("chr") else c, p, r, a)


def spliceai_ds(info_value):
    """{gene (version-stripped): DS_MAX} from a SpliceAI v1 INFO value: 9-field '|' blocks, one per gene."""
    s = info_value.split("|"); out = {}
    for b in range(len(s) // 9):
        gf = s[b * 9 + 1].split("---")
        gid = (gf[1] if len(gf) > 1 else gf[0]).split(".")[0]
        try:
            ds = max(float(x) for x in s[b * 9 + 2:b * 9 + 6])
        except ValueError:
            continue
        out[gid] = max(ds, out.get(gid, -1.0))
    return out


def read_spliceai_v1(path):
    ds = {}
    with opn(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t", 8)
            vid = norm("%s:%s:%s:%s" % (f[0], f[1], f[3], f[4]))
            for kv in f[7].split(";"):
                if kv.startswith("SpliceAI="):
                    for g, v in spliceai_ds(kv[9:]).items():
                        ds[(vid, g)] = max(v, ds.get((vid, g), -1.0))
    return ds


def read_spliceai2(paths):
    """{(variant key, strand): score} from SpliceAI2 tables: chrom pos ref alt strand ... score (last column);
    a header line, if any, is skipped. Several files: precomputed lookups and GPU scores."""
    s2 = {}
    for p in paths:
        with opn(p) as fh:
            for l in fh:
                f = l.rstrip("\n").split("\t")
                if len(f) < 6 or not f[1].isdigit():
                    continue
                s2[(norm("%s:%s:%s:%s" % tuple(f[:4])), f[4])] = float(f[-1])
    return s2


def read_gtf(path):
    """GENCODE gene strand and span: {gene: strand}, {(gene, chrom without 'chr'): (start, end)}."""
    strand, span = {}, {}
    with opn(path) as fh:
        for l in fh:
            if l[0] == "#":
                continue
            f = l.split("\t", 8)
            if f[2] != "gene":
                continue
            g = re.search(r'gene_id "([^"]+)"', f[8]).group(1).split(".")[0]
            strand[g] = f[6]
            c = f[0][3:] if f[0].startswith("chr") else f[0]
            span[(g, c)] = (int(f[3]), int(f[4]))   # per chromosome: PAR genes sit on X and Y
    return strand, span


def read_cadd_indels(path):
    """CADD's own output (Chrom Pos Ref Alt RawScore PHRED) -> {variant key: PHRED}."""
    out = {}
    with opn(path) as fh:
        for l in fh:
            if l[0] == "#":
                continue
            f = l.rstrip("\n").split("\t")
            out[norm("%s:%s:%s:%s" % tuple(f[:4]))] = float(f[5])
    return out


def fmt(x):
    return "NA" if x != x else "%g" % x


def run(a, log=sys.stdout):
    model = "spliceai2" if a.spliceai2 else "spliceai"
    T = MODELS[model]
    n = Counter()
    if model == "spliceai2":
        if not a.gtf:
            sys.exit("--spliceai2 needs --gtf (GENCODE v39): an entry takes the score on its gene's strand, inside its gene")
        s2 = read_spliceai2(a.spliceai2); strand, span = read_gtf(a.gtf)
    else:
        v1 = read_spliceai_v1(a.spliceai)
    cadd_ind = read_cadd_indels(a.cadd_indels) if a.cadd_indels else {}

    # One row per (variant, gene): the MANE Select row, else the canonical one, protein-coding only.
    keep = OrderedDict()
    with opn(a.vep) as fh:
        h = fh.readline().split()
        need = [a.vep_snp_id_col, a.vep_gene_col, a.vep_lof_col, a.vep_revel_col, a.vep_cadd_phred_col,
                a.vep_consequence_col, a.vep_transcript_col, a.vep_mane_select_col, a.vep_canonical_col, a.vep_biotype_col]
        miss = [c for c in need if c not in h]
        if miss:
            sys.exit("VEP table lacks column(s) %s (header: %s)" % (miss, " ".join(h)))
        ix = [h.index(c) for c in need]
        for line in fh:
            f = line.split()
            if not f:
                continue
            vid, gene, lof, revel, cadd, csq, tx, mane, canon, bio = (f[i] if i < len(f) else "." for i in ix)
            n["VEP rows"] += 1
            mane = "" if mane in (".", "NA", "") else mane
            if bio != "protein_coding" or not (mane or canon == "YES"):
                continue
            k = (vid, gene.split(".")[0])
            row = dict(LOF=lof, REVEL=revel, CADD=cadd, CSQ=csq, TRANSCRIPT=tx, MANE=mane, CANON=canon)
            if k in keep:
                n["(variant, gene) with two MANE/canonical rows: MANE Select row kept"] += 1
                if keep[k]["MANE"] or not mane:
                    continue
            keep[k] = row
    n["(variant, gene) entries"] = len(keep)

    out_long = a.out_file + ".long.tsv.gz"
    groups = OrderedDict(); labels = Counter(); src = Counter(); lofmiss = []
    with gzip.open(out_long, "wt") as fl:
        fl.write("\t".join(LONG_COLS) + "\n")
        for (vid, gene), r in keep.items():
            key = norm(vid)
            lof = r["LOF"] if r["LOF"] in ("HC", "LC") else ""
            if not lof and any(c in LOFTEE_CSQS for c in r["CSQ"].split("&")):
                lofmiss.append(vid)
            revel, cadd = first_value(r["REVEL"]), first_value(r["CADD"])
            if cadd != cadd and key in cadd_ind:
                cadd = cadd_ind[key]; n["CADD from --cadd_indels"] += 1
            if model == "spliceai2":
                v = s2.get((key, strand.get(gene, "?")))
                if v is None:
                    ds, sc = NAN, "none"
                else:
                    s_, e_ = span.get((gene, key.split(":")[0]), (0, -1))
                    sc = "spliceai2" if s_ <= int(key.split(":")[1]) <= e_ else "spliceai2_outside_gene"
                    ds = v
            else:
                ds = v1.get((key, gene), NAN); sc = "spliceai" if ds == ds else "none"
            src[sc] += 1
            lab = annotate(r["CSQ"], lof, revel, cadd, ds if sc in ("spliceai", "spliceai2") else NAN,
                           T["t_hi"], T["t_lo"], T["label"])
            labels[lab] += 1
            groups.setdefault(gene, []).append((vid, lab))
            fl.write("\t".join([vid, gene, lof or "NA", fmt(revel), fmt(cadd), r["CSQ"], r["TRANSCRIPT"],
                                r["MANE"] or "NA", r["CANON"] or "NA", "protein_coding", fmt(ds), sc, lab]) + "\n")
    if lofmiss and not a.allow_missing_loftee:
        os.remove(out_long)
        sys.exit("%d entries have a stop_gained / frameshift / splice donor or acceptor consequence but no LOFTEE call "
                 "(e.g. %s). LOFTEE calls every such entry HC or LC, so it failed on these variants: check VEP's "
                 "warnings file for 'Plugin 'LoF' went wrong', re-run VEP on them without --fork, and merge. "
                 "--allow_missing_loftee writes the files anyway." % (len(lofmiss), ", ".join(lofmiss[:3])))
    if keep and src["spliceai"] + src["spliceai2"] + src["spliceai2_outside_gene"] == 0:
        os.remove(out_long)
        sys.exit("No SpliceAI score matched any entry: check that the variant IDs (chrN:pos:ref:alt), the chromosome "
                 "names and, for SpliceAI2, the GTF all describe the same build and the same variants.")
    with opn(a.out_file, "wt") as fo:
        for gene, vs in groups.items():
            fo.write("%s var %s\n%s anno %s\n" % (gene, " ".join(v for v, _ in vs), gene, " ".join(l for _, l in vs)))
    for k, v in n.items():
        print("  %10d  %s" % (v, k), file=log)
    for k, v in src.most_common():
        print("  %10d  SpliceAI source: %s" % (v, k), file=log)
    for k in ["pLoF", "damaging_missense", T["label"], "other_missense", "synonymous", "non_coding"]:
        print("  %10d  %s" % (labels[k], k), file=log)
    print("wrote %s (%d genes) and %s" % (a.out_file, len(groups), out_long), file=log)
    return labels


# ------------------------------------------------------------------------------------------------ selftest
def selftest():
    ok = True

    def check(name, cond):
        nonlocal ok
        print("  %-76s %s" % (name, "PASS" if cond else "FAIL")); ok &= bool(cond)

    S2 = MODELS["spliceai2"]; V1 = MODELS["spliceai"]
    A2 = lambda csq, lof="", revel=NAN, cadd=NAN, ds=NAN: annotate(csq, lof, revel, cadd, ds, S2["t_hi"], S2["t_lo"], S2["label"])
    A1 = lambda csq, lof="", revel=NAN, cadd=NAN, ds=NAN: annotate(csq, lof, revel, cadd, ds, V1["t_hi"], V1["t_lo"], V1["label"])
    check("HC -> pLoF, whatever else", A2("stop_gained", "HC", ds=0.9) == "pLoF")
    check("LC -> damaging_missense", A2("frameshift_variant", "LC") == "damaging_missense")
    check("missense REVEL 0.773 -> damaging; 0.772 and CADD 28.0 -> other_missense",
          A2("missense_variant", revel=0.773) == "damaging_missense" and A2("missense_variant", revel=0.772, cadd=28.0) == "other_missense")
    check("protein_altering_variant CADD 28.1 -> damaging", A2("protein_altering_variant", cadd=28.1) == "damaging_missense")
    check("damaging before the splice category", A2("missense_variant", revel=0.9, ds=0.9) == "damaging_missense")
    check("SpliceAI2 0.40 -> spliceai2_gt040; 0.39 intron -> non_coding",
          A2("intron_variant", ds=0.40) == "spliceai2_gt040" and A2("intron_variant", ds=0.39) == "non_coding")
    check("SpliceAI v1 0.5 -> spliceai_gt05; 0.49 -> non_coding",
          A1("intron_variant", ds=0.5) == "spliceai_gt05" and A1("intron_variant", ds=0.49) == "non_coding")
    check("missense low REVEL with SpliceAI2 0.5 -> splice category",
          A2("missense_variant", revel=0.1, cadd=5, ds=0.5) == "spliceai2_gt040")
    check("missense SpliceAI2 0.3 -> other_missense (no route into damaging)", A2("missense_variant", ds=0.3) == "other_missense")
    check("synonymous: < 0.10 synonymous, 0.10-0.40 non_coding, >= 0.40 splice",
          A2("synonymous_variant", ds=0.09) == "synonymous" and A2("synonymous_variant", ds=0.10) == "non_coding"
          and A2("synonymous_variant", ds=0.40) == "spliceai2_gt040")
    check("synonymous without a score -> synonymous", A2("synonymous_variant") == "synonymous")
    check("splice_region_variant&synonymous_variant -> synonymous (any term)",
          A2("splice_region_variant&synonymous_variant") == "synonymous")
    check("synonymous_variant&intron_variant -> non_coding (non-coding tested first)",
          A2("synonymous_variant&intron_variant") == "non_coding")
    check("missense_variant&splice_region_variant counts as missense",
          A2("missense_variant&splice_region_variant", revel=0.8) == "damaging_missense")
    check("unmatched consequence (splice_region alone) -> non_coding", A2("splice_region_variant") == "non_coding")
    check("REVEL lists: '&', ';', ',', '%3B' and '.'", first_value(".&0.81&0.2") == 0.81 and
          first_value(".;0.53") == first_value(".,0.53") == first_value(".%3B0.53") == 0.53 and first_value(".") != first_value("."))
    sp = ("T|RPL36---ENSG00000130255.13---ENST1---yes---protein_coding---NM_1,NM_2|0.10|0.20|0.30|0.40|-6|-27|6|-27,"
          "T|GB---ENSG00000999999.1---ENST2---yes---protein_coding---NM_3|0.90|0.00|0.00|0.00|21|2|1|8")
    check("SpliceAI v1 parse survives a RefSeq comma list", spliceai_ds(sp) == {"ENSG00000130255": 0.4, "ENSG00000999999": 0.9})
    check("variant keys join across chr prefix and separators", norm("chr1:100:A:G") == norm("1_100_A_G") == "1:100:A:G")

    # end to end, SpliceAI2 then v1, on a tiny VEP table
    t = tempfile.mkdtemp(); null = open(os.devnull, "w")
    vep = os.path.join(t, "vep.txt")
    rows_in = [
        "chr1:100:A:G ENSG1 . .&0.9 20 missense_variant ENST1 NM_1 YES protein_coding",
        "chr1:100:A:G ENSG1 . . . missense_variant ENST1b . . protein_coding",          # not MANE/canonical
        "chr1:100:A:G ENSG2 . . . upstream_gene_variant ENST2 . YES lncRNA",          # not protein-coding
        "chr1:200:C:T ENSG1 . . . intron_variant ENST1 NM_1 YES protein_coding",       # SpliceAI2 0.6 -> splice
        "chr1:300:G:A ENSG1 . . . synonymous_variant ENST1 NM_1 YES protein_coding",   # 0.2 -> non_coding
        "chr1:400:T:C ENSG1 HC . . stop_gained ENST1 NM_1 YES protein_coding",
        "chr1:900:T:C ENSG1 . . . downstream_gene_variant ENST1 NM_1 YES protein_coding",  # outside ENSG1: no route
        "chr1:950:G:T ENSG1 . . . synonymous_variant ENST1 . YES protein_coding",      # canonical, then MANE row wins
        "chr1:950:G:T ENSG1 . . . missense_variant ENST1m NM_1 . protein_coding",
    ]
    with open(vep, "w") as fh:
        fh.write("SNP_ID GENE LOF REVEL_SCORE CADD_PHRED CSQ TRANSCRIPT MANE_SELECT CANONICAL BIOTYPE\n")
        fh.write("\n".join(rows_in) + "\n")
    gtf = os.path.join(t, "g.gtf")
    with open(gtf, "w") as fh:
        fh.write('chr1\tHAVANA\tgene\t50\t800\t.\t+\t.\tgene_id "ENSG1.3";\n')
    s2 = os.path.join(t, "s2.tsv")
    with open(s2, "w") as fh:
        fh.write("chrom\tpos\tref\talt\tstrand\tspliceai2_summary_score\n")
        for r in ["chr1\t200\tC\tT\t+\t0.6", "chr1\t200\tC\tT\t-\t0.0", "chr1\t300\tG\tA\t+\t0.2",
                  "chr1\t900\tT\tC\t+\t0.9", "chr1\t100\tA\tG\t-\t0.99"]:   # 100 '-': the other strand, not ENSG1's
            fh.write(r + "\n")
    out = os.path.join(t, "o.txt.gz")
    ns = argparse.Namespace(vep=vep, spliceai=None, spliceai2=[s2], gtf=gtf, out_file=out, cadd_indels=None,
                            allow_missing_loftee=False, **DEFAULT_COLS)
    run(ns, log=null)
    rows = {tuple(l.split("\t")[:2]): l.rstrip("\n").split("\t") for l in list(gzip.open(out + ".long.tsv.gz", "rt"))[1:]}
    lab = {k: v[-1] for k, v in rows.items()}
    check("end to end: MANE/canonical protein-coding rows only (6 entries)", len(rows) == 6)
    check("end to end: REVEL '.&0.9' -> damaging_missense", lab[("chr1:100:A:G", "ENSG1")] == "damaging_missense")
    check("end to end: SpliceAI2 score taken on the gene's strand only", rows[("chr1:100:A:G", "ENSG1")][10] == "NA")
    check("end to end: intron 0.6 -> spliceai2_gt040", lab[("chr1:200:C:T", "ENSG1")] == "spliceai2_gt040")
    check("end to end: synonymous 0.2 -> non_coding", lab[("chr1:300:G:A", "ENSG1")] == "non_coding")
    check("end to end: HC -> pLoF", lab[("chr1:400:T:C", "ENSG1")] == "pLoF")
    check("end to end: outside the gene -> no route, score kept",
          lab[("chr1:900:T:C", "ENSG1")] == "non_coding" and rows[("chr1:900:T:C", "ENSG1")][10:12] == ["0.9", "spliceai2_outside_gene"])
    check("end to end: two rows for one (variant, gene) -> the MANE Select row",
          rows[("chr1:950:G:T", "ENSG1")][6] == "ENST1m" and lab[("chr1:950:G:T", "ENSG1")] == "other_missense")
    g = [l.split() for l in gzip.open(out, "rt")]
    check("end to end: group file has var/anno lines in long-table order",
          g[0][:2] == ["ENSG1", "var"] and g[1][:2] == ["ENSG1", "anno"] and len(g[0]) == len(g[1]) == 8
          and g[1][2:] == [lab[(v, "ENSG1")] for v in g[0][2:]])
    v1 = os.path.join(t, "v1.vcf")
    with open(v1, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        fh.write("1\t200\t.\tC\tT\t.\t.\tSpliceAI=T|X---ENSG1.5---ENST1---yes---protein_coding---NM_1|0.55|0.00|0.00|0.00|1|2|3|4\n")
        fh.write("1\t300\t.\tG\tA\t.\t.\tSpliceAI=A|X---ENSG1.5---ENST1---yes---protein_coding---NM_1|0.30|0.00|0.00|0.00|1|2|3|4\n")
    ns.spliceai, ns.spliceai2, ns.gtf, ns.out_file = v1, None, None, os.path.join(t, "o1.txt")
    run(ns, log=null)
    lab1 = {tuple(l.split("\t")[:2]): l.rstrip("\n").split("\t")[-1] for l in list(gzip.open(ns.out_file + ".long.tsv.gz", "rt"))[1:]}
    check("end to end, SpliceAI v1 (no 'chr' in its VCF): 0.55 -> spliceai_gt05, synonymous 0.3 -> non_coding",
          lab1[("chr1:200:C:T", "ENSG1")] == "spliceai_gt05" and lab1[("chr1:300:G:A", "ENSG1")] == "non_coding")

    def refused(ns_):
        try:
            run(ns_, log=null); return False
        except SystemExit:
            return True
    nn = argparse.Namespace(**vars(ns)); nn.spliceai = os.path.join(t, "none.vcf"); nn.out_file = os.path.join(t, "o2.txt")
    open(nn.spliceai, "w").write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    check("negative control: no SpliceAI score matching any entry is refused", refused(nn))
    with open(vep, "a") as fh:   # a stop_gained with no LOFTEE call: LOFTEE failed on it
        fh.write("chr1:500:C:A ENSG1 . . . stop_gained ENST1 NM_1 YES protein_coding\n")
    ns.out_file = os.path.join(t, "o3.txt")
    check("negative control: a stop_gained without a LOFTEE call is refused", refused(ns))
    ns.allow_missing_loftee = True
    check("... and written with --allow_missing_loftee", not refused(ns))
    print("\nselftest %s" % ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vep", "-v", required=True, help="VEP table (step 1; whitespace-delimited, header line)")
    sp = ap.add_mutually_exclusive_group(required=True)
    sp.add_argument("--spliceai2", nargs="+", help="SpliceAI2 score table(s): chrom pos ref alt strand ... score (default model)")
    sp.add_argument("--spliceai", "-s", help="SpliceAI v1 VCF instead of SpliceAI2 (label spliceai_gt05)")
    ap.add_argument("--gtf", help="GENCODE v39 GTF (with --spliceai2): gene strand and span")
    ap.add_argument("--out_file", "-w", required=True, help="group file; the long table goes to OUT.long.tsv.gz")
    ap.add_argument("--cadd_indels", help="CADD output for indels (step 3), used where VEP has no CADD")
    ap.add_argument("--allow_missing_loftee", action="store_true",
                    help="write the files even when LoF-type entries lack a LOFTEE call (otherwise refused)")
    for k, v in DEFAULT_COLS.items():
        ap.add_argument("--" + k, default=v, help="column name in the VEP table (default %s)" % v)
    run(ap.parse_args())

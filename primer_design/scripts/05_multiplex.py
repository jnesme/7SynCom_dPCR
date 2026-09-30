#!/usr/bin/env python
"""Choose one set per strain that is mutually compatible for multiplexing.

Pool = best set per gene among specificity-passing sets (mid-replichore first, then
route class, fewest 3-mismatch background sites, primer3 penalty); top N genes per strain.
Pairwise interactions between sets (primer3 thermodynamics, design buffer, 37 C):
  any-dimer dG between every oligo pair of two different sets (primers + probes),
  3'-anchored dG where the extendable 3' end is a primer (probes are 3'-blocked).
Within-set dimers were constrained by primer3 at design time and are not part of the
between-set objective. All 7-way combinations are enumerated; objective (lexicographic):
  max worst any-dimer dG -> max worst 3'-dimer dG -> min primer Tm spread -> min class_rank sum.
Best combos are checked for cross-set products (exhaustive in-silico PCR with every primer
of the pool, insilico_pcr.py) and for probes binding other sets' amplicons.
Also proposes a 4+3 split into two wells for <=5-colour instruments.
"""
import itertools, os
import numpy as np
import pandas as pd
import primer3
from insilico_pcr import build_index, find_sites, products

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
RES = os.path.join(HERE, "..", "results")
THREADS = os.environ.get("THREADS", "1")
N_PER_STRAIN = 6
MAX_PROD = 3000
PROBE_MIN_MM = 5  # probe must have >= 5 mismatches to any other amplicon
COND = dict(mv_conc=50.0, dv_conc=3.8, dntp_conc=0.8, dna_conc=800.0, temp_c=37.0)
allg = os.path.join(WORK, "all_genomes.fna")
idx = os.path.join(WORK, "bt_all")
rc = str.maketrans("ACGT", "TGCA")

rep = pd.read_csv(os.path.join(RES, "specificity_report.tsv"), sep="\t")
ok = rep[rep["pass"]]
ok = ok.assign(mid=ok.rel_ori.between(0.3, 0.7)).sort_values(
    ["strain", "mid", "class_rank", "n_3mm_sites", "penalty"], ascending=[True, False, True, True, True])
missing = sorted(set(rep.strain) - set(ok.strain))
if missing:
    raise SystemExit(f"no specificity-passing set for {missing}; relax design or add fallback loci")
pool = ok.groupby(["strain", "locus_tag"], sort=False).head(1).groupby("strain").head(N_PER_STRAIN).reset_index(drop=True)
strains = sorted(pool.strain.unique())
print(pool.groupby("strain").agg(n=("set_id", "size"), mid_replichore=("mid", "sum")))

oligos = {i: [("fwd", r.fwd), ("rev", r.rev), ("probe", r.probe)] for i, r in pool.iterrows()}
dg_cache, end_cache = {}, {}


def dg_any(a, b):
    k = (a, b) if a < b else (b, a)
    if k not in dg_cache:
        dg_cache[k] = primer3.calc_heterodimer(a, b, **COND).dg / 1000.0
    return dg_cache[k]


def dg_end(a, b):
    """3' end of a annealed on b (asymmetric)."""
    if (a, b) not in end_cache:
        end_cache[(a, b)] = primer3.calc_end_stability(a, b, **COND).dg / 1000.0
    return end_cache[(a, b)]


n = len(pool)
ANY = np.zeros((n, n)); END = np.zeros((n, n))  # diagonal (within-set) stays 0 = neutral
for i, j in itertools.combinations(range(n), 2):
    if pool.strain[i] == pool.strain[j]:
        continue
    pairs = [(x, y) for x in oligos[i] for y in oligos[j]]
    ANY[i, j] = ANY[j, i] = min(dg_any(x[1], y[1]) for x, y in pairs)
    ends = [dg_end(x[1], y[1]) for x, y in pairs if x[0] != "probe"] + \
           [dg_end(y[1], x[1]) for x, y in pairs if y[0] != "probe"]
    END[i, j] = END[j, i] = min(ends)

tms = pool[["fwd_tm", "rev_tm"]].values
by_strain = [list(pool.index[pool.strain == s]) for s in strains]
res = []
for combo in itertools.product(*by_strain):
    c = np.array(combo)
    t = tms[c].ravel()
    res.append((ANY[np.ix_(c, c)].min(), END[np.ix_(c, c)].min(), t.max() - t.min(),
                pool.class_rank[c].sum(), combo))
res.sort(key=lambda x: (-round(x[0], 1), -round(x[1], 1), round(x[2], 1), x[3]))
print(f"{len(res):,} combinations scored; best worst-dG any={res[0][0]:.2f} end={res[0][1]:.2f} kcal/mol")

# binding sites of every pool primer, computed once (<= 3 mismatches, exhaustive)
build_index(allg, idx, THREADS)
psites = find_sites({f"{i}__{role}": seq for i in pool.index for role, seq in oligos[i] if role != "probe"},
                    idx, WORK, "pool", threads=THREADS)
psites["set"] = psites.name.str.split("__").str[0].astype(int)


def cross_check(combo):
    sel = pool.loc[list(combo)]
    p = products(psites[psites.set.isin(combo)], max_len=MAX_PROD)
    p = p.drop_duplicates(["chrom", "start", "end"])
    expected = {(f"{r.strain}|{r.seqid}", r.amp_start - 1, r.amp_end) for r in sel.itertuples()}
    extra = p[[k not in expected for k in zip(p.chrom, p.start, p.end)]]
    found = set(zip(p.chrom, p.start, p.end))
    probe_hits = []  # exhaustive ungapped scan of each probe on both strands of the other amplicons
    for a in sel.itertuples():
        L = len(a.probe)
        for b in sel.itertuples():
            if a.set_id == b.set_id:
                continue
            for t in (b.amplicon, b.amplicon.translate(rc)[::-1]):
                mm = min(sum(x != y for x, y in zip(a.probe, t[i:i + L])) for i in range(len(t) - L + 1))
                if mm < PROBE_MIN_MM:
                    probe_hits.append((a.set_id, b.set_id, mm))
    return expected <= found, extra, probe_hits


chosen = None
for k, r in enumerate(res[:500]):
    ok_on, extra, probe_hits = cross_check(r[4])
    if ok_on and extra.empty and not probe_hits:
        chosen = r; break
    print(f"combo #{k} rejected: on-target={ok_on}, extra products={len(extra)}, probe hits={len(probe_hits)}")
if chosen is None:
    raise SystemExit("no combination passed the cross-reactivity check")

sel = pool.loc[list(chosen[4])].copy()
print(f"\nchosen combination: worst between-set any-dimer dG {chosen[0]:.2f}, worst 3'-dimer dG "
      f"{chosen[1]:.2f} kcal/mol, primer Tm spread {chosen[2]:.2f} C")

# 4+3 well split maximising the worst within-well (between-set) interaction
ids = list(sel.index)
best = None
for w1 in itertools.combinations(ids, 4):
    w2 = [i for i in ids if i not in w1]
    score = min(ANY[np.ix_(w, w)].min() for w in (list(w1), w2))
    if best is None or score > best[0]:
        best = (score, w1, w2)
sel["well_4plus3"] = ["W1" if i in best[1] else "W2" for i in sel.index]

strain_info = pd.read_csv(os.path.join(WORK, "strains.tsv"), sep="\t").set_index("strain")
sel.insert(1, "species_NCBI", sel.strain.map(strain_info.NCBI_species))
sel.to_csv(os.path.join(RES, "final_multiplex.tsv"), sep="\t", index=False)
pool.to_csv(os.path.join(RES, "candidates_pool.tsv"), sep="\t", index=False)
ok.to_csv(os.path.join(RES, "candidates_all.tsv"), sep="\t", index=False)

ol = [(f"{r.strain}_{role}", getattr(r, role)) for r in sel.itertuples() for role in ("fwd", "rev", "probe")]
m = pd.DataFrame([[dg_any(a, b) for _, b in ol] for _, a in ol], index=[x for x, _ in ol], columns=[x for x, _ in ol])
m.round(2).to_csv(os.path.join(RES, "dimer_matrix.tsv"), sep="\t")
print(sel[["strain", "species_NCBI", "locus_tag", "gene", "product", "rel_ori", "amp_len", "fwd", "rev", "probe",
           "well_4plus3"]].to_string(index=False))

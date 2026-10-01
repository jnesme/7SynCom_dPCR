#!/usr/bin/env python
"""ADD-ON: multiplex selection that minimises the strongest dimer of ANY kind.

Does not replace step 05. Step 05 takes the 6 best-ranked genes per strain and only
optimises dimers BETWEEN assays, so (a) within-assay dimers and self-dimers were left
to primer3's design thresholds and (b) the worst between-assay pair is capped by that
small pool. This add-on instead:

  1. keeps only assays on strain-unique genes (homology) that passed specificity;
  2. scores every assay on its own oligos: self-dimers F-F, R-R, P-P and F-R, F-P, R-P;
  3. takes up to N_PER_STRAIN genes per strain (best assay per gene) - a much wider pool;
  4. finds the 7-plex (one assay per strain) whose WORST oligo pair, within or between
     assays, is as weak as possible (max-min heterodimer dG, primer3, design buffer,
     800 nM, 37 C). Exact search: binary search on the dG threshold + backtracking.
     Ties are resolved toward better-ranked assays (mid-replichore, route class, fewer
     3-mismatch background sites);
  5. checks the winner like step 05 (pooled exhaustive in-silico PCR, probes vs the
     other amplicons) and proposes a 4+3 well split.

Outputs go to primer_design/strict/results/, so the standard results are untouched.
run_strict_addon.sh then runs the unmodified steps 06 and 07 on that folder.
"""
import itertools, os
import numpy as np
import pandas as pd
import primer3
from insilico_pcr import build_index, find_sites, products

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
BASE_RES = os.path.join(HERE, "..", "results")          # standard run (input)
RES = os.path.join(HERE, "..", "strict", "results")      # add-on output
os.makedirs(RES, exist_ok=True)
THREADS = os.environ.get("THREADS", "1")
N_PER_STRAIN = 40      # genes per strain in the pool (step 05 uses 6)
MAX_PROD = 3000
PROBE_MIN_MM = 5       # probe must have >= 5 mismatches to any other amplicon
COND = dict(mv_conc=50.0, dv_conc=3.8, dntp_conc=0.8, dna_conc=800.0, temp_c=37.0)
allg = os.path.join(WORK, "all_genomes.fna")
idx = os.path.join(WORK, "bt_all")
rc = str.maketrans("ACGT", "TGCA")

rep = pd.read_csv(os.path.join(BASE_RES, "specificity_report.tsv"), sep="\t")
rep.to_csv(os.path.join(RES, "specificity_report.tsv"), sep="\t", index=False)  # read by steps 06/07
ok = rep[rep["pass"] & rep.hom_unique].copy()
missing = sorted(set(rep.strain) - set(ok.strain))
if missing:
    raise SystemExit(f"no specificity-passing assay on a strain-unique gene for {missing}")

dg_cache = {}


def dimer(a, b):
    """(dG kcal/mol, Tm) of the most stable duplex between two oligos."""
    k = (a, b) if a < b else (b, a)
    if k not in dg_cache:
        t = primer3.calc_heterodimer(a, b, **COND)
        dg_cache[k] = (t.dg / 1000.0, t.tm)
    return dg_cache[k]


def within_assay(r):
    """Strongest dimer among the assay's own oligos, self-dimers included."""
    d = [dimer(a, b) for a, b in itertools.combinations_with_replacement([r.fwd, r.rev, r.probe], 2)]
    return pd.Series({"within_dg": min(x[0] for x in d), "within_tm": max(x[1] for x in d)})


ok[["within_dg", "within_tm"]] = ok.apply(within_assay, axis=1)
ok = ok.assign(mid=ok.rel_ori.between(0.3, 0.7)).sort_values(
    ["strain", "mid", "class_rank", "n_3mm_sites", "penalty"], ascending=[True, False, True, True, True])
# best assay per gene = the one with the weakest within-assay dimer
best = ok.sort_values("within_dg", ascending=False).groupby(["strain", "locus_tag"], sort=False).head(1)
pool = (best.sort_values(["strain", "mid", "class_rank", "n_3mm_sites", "penalty"],
                         ascending=[True, False, True, True, True])
        .groupby("strain").head(N_PER_STRAIN).reset_index(drop=True))
strains = sorted(pool.strain.unique())
print(pool.groupby("strain").agg(genes=("set_id", "size"), mid_replichore=("mid", "sum"),
                                 best_within_dg=("within_dg", "max")).round(2))

n = len(pool)
olig = [[r.fwd, r.rev, r.probe] for r in pool.itertuples()]
BETW = np.zeros((n, n))
for i, j in itertools.combinations(range(n), 2):
    if pool.strain[i] != pool.strain[j]:
        BETW[i, j] = BETW[j, i] = min(dimer(x, y)[0] for x in olig[i] for y in olig[j])
within = pool.within_dg.values
by_strain = [list(pool.index[pool.strain == s]) for s in strains]
# search the most constrained strain first
order = sorted(range(len(strains)), key=lambda k: len(by_strain[k]))


def solutions(t):
    """Yield 7-plexes whose every oligo pair is weaker than t (dG > t), best-ranked first."""
    cand = [[i for i in by_strain[k] if within[i] > t] for k in order]
    if any(len(c) == 0 for c in cand):
        return

    def rec(k, chosen):
        if k == len(cand):
            yield tuple(chosen); return
        for i in cand[k]:
            if all(BETW[i, j] > t for j in chosen):
                yield from rec(k + 1, chosen + [i])
    yield from rec(0, [])


def cross_check(combo):
    sel = pool.loc[list(combo)]
    p = products(psites[psites.set.isin(combo)], max_len=MAX_PROD).drop_duplicates(["chrom", "start", "end"])
    expected = {(f"{r.strain}|{r.seqid}", r.amp_start - 1, r.amp_end) for r in sel.itertuples()}
    found = set(zip(p.chrom, p.start, p.end))
    extra = [k for k in found if k not in expected]
    probe_hits = []
    for a in sel.itertuples():
        L = len(a.probe)
        for b in sel.itertuples():
            if a.set_id == b.set_id:
                continue
            for t in (b.amplicon, b.amplicon.translate(rc)[::-1]):
                mm = min(sum(x != y for x, y in zip(a.probe, t[i:i + L])) for i in range(len(t) - L + 1))
                if mm < PROBE_MIN_MM:
                    probe_hits.append((a.set_id, b.set_id, mm))
    return expected <= found and not extra and not probe_hits


build_index(allg, idx, THREADS)
psites = find_sites({f"{i}__{role}": seq for i in pool.index for role, seq in zip(("fwd", "rev"), olig[i][:2])},
                    idx, WORK, "pool_strict", threads=THREADS)
psites["set"] = psites.name.str.split("__").str[0].astype(int)

# binary search on the threshold over the distinct dG values
levels = np.unique(np.round(np.concatenate([within, BETW[BETW < 0]]), 2))
lo, hi, chosen, chosen_t = 0, len(levels) - 1, None, None
while lo <= hi:
    mid_i = (lo + hi) // 2
    t = levels[mid_i] - 0.005
    sol = next((c for c in itertools.islice(solutions(t), 2000) if cross_check(c)), None)
    if sol is not None:
        chosen, chosen_t = sol, t
        lo = mid_i + 1      # try a stricter (less negative) threshold
    else:
        hi = mid_i - 1
if chosen is None:
    raise SystemExit("no combination found")

sel = pool.loc[sorted(chosen, key=lambda i: pool.strain[i])].copy()
c = list(sel.index)
worst_between = min(BETW[i, j] for i, j in itertools.combinations(c, 2))
print(f"\nchosen 7-plex: worst pair of any kind {min(worst_between, sel.within_dg.min()):.2f} kcal/mol "
      f"(between assays {worst_between:.2f}, within an assay {sel.within_dg.min():.2f}); "
      f"highest dimer Tm within an assay {sel.within_tm.max():.1f} C")

ids = list(sel.index)
bestsplit = None
for w1 in itertools.combinations(ids, 4):
    w2 = [i for i in ids if i not in w1]
    score = min(min(BETW[i, j] for i, j in itertools.combinations(w, 2)) for w in (list(w1), w2))
    if bestsplit is None or score > bestsplit[0]:
        bestsplit = (score, w1, w2)
sel["well_4plus3"] = ["W1" if i in bestsplit[1] else "W2" for i in sel.index]

strain_info = pd.read_csv(os.path.join(WORK, "strains.tsv"), sep="\t").set_index("strain")
sel.insert(1, "species_NCBI", sel.strain.map(strain_info.NCBI_species))
sel.to_csv(os.path.join(RES, "final_multiplex.tsv"), sep="\t", index=False)
# step 07 (figure 4B) enumerates all combinations of candidates_pool.tsv: give it the
# 6 best-ranked genes per strain of this pool, always including the chosen assays
small = pd.concat([sel.drop(columns=["species_NCBI", "well_4plus3"]),
                   pool.drop(index=sel.index).groupby("strain").head(5)]).sort_values("strain")
small.to_csv(os.path.join(RES, "candidates_pool.tsv"), sep="\t", index=False)
pool.to_csv(os.path.join(RES, "candidates_pool_wide.tsv"), sep="\t", index=False)
# worst between-assay dG of every combination of that reduced pool, chosen set first
# (read by the R figure script, see R/run_addon_figures_R.sh)
grp = [list(small.index[small.strain == s]) for s in strains]
sc = pd.DataFrame([(min(BETW[i, j] for i, j in itertools.combinations(cb, 2)), set(cb) == set(c))
                   for cb in itertools.product(*grp)], columns=["worst_any", "chosen"])
sc.sort_values(["chosen", "worst_any"], ascending=False).to_csv(
    os.path.join(RES, "combination_scores.tsv"), sep="\t", index=False)
ok.to_csv(os.path.join(RES, "candidates_all.tsv"), sep="\t", index=False)

ol = [(f"{r.strain}_{role}", getattr(r, role)) for r in sel.itertuples() for role in ("fwd", "rev", "probe")]
m = pd.DataFrame([[dimer(a, b)[0] for _, b in ol] for _, a in ol], index=[x for x, _ in ol], columns=[x for x, _ in ol])
m.round(2).to_csv(os.path.join(RES, "dimer_matrix.tsv"), sep="\t")
print(sel[["strain", "species_NCBI", "locus_tag", "gene", "product", "rel_ori", "amp_len", "within_dg",
           "n_3mm_sites", "fwd", "rev", "probe", "well_4plus3"]].to_string(index=False))

#!/usr/bin/env python
"""Publication figures (PNG + PDF) in primer_design/figures/.

Fig 1  workflow with the number of genes / assays surviving each step
Fig 2  primary chromosome maps: unique sequence, passing assays, chosen target, ori/ter
Fig 3  in-silico specificity: assay x genome matrix + knock-out test
Fig 4  multiplex compatibility: oligo heterodimer matrix + score of all 7-way combinations

Colour: sequential blue ramp for magnitudes, one orange accent for the chosen assays,
neutral greys otherwise (palette validated with the dataviz validator, light mode).
"""
import itertools, os
import numpy as np
import pandas as pd
import primer3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle, Wedge
from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)

# ---- tokens --------------------------------------------------------------------
SURF, INK, INK2, MUTED, GRID = "#ffffff", "#0b0b0b", "#52514e", "#8a8983", "#e4e3df"
BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#b9b8b3"
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("seq_blue", RAMP)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "figure.facecolor": SURF,
    "axes.facecolor": SURF, "savefig.facecolor": SURF, "pdf.fonttype": 42, "svg.fonttype": "none",
})
COND = dict(mv_conc=50.0, dv_conc=3.8, dntp_conc=0.8, dna_conc=800.0, temp_c=37.0)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def panel_letter(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="right")


# ---- shared data ------------------------------------------------------------------
strains = pd.read_csv(os.path.join(WORK, "strains.tsv"), sep="\t").set_index("strain")
species = strains.NCBI_species.to_dict()
cds = pd.read_csv(os.path.join(WORK, "cds_classified.tsv"), sep="\t")
txt = pd.read_csv(os.path.join(WORK, "cds_text.tsv"), sep="\t").fillna("")
cds = cds.merge(txt, on="id")
cds["class"] = cds["class"].fillna("")
cand = pd.read_csv(os.path.join(WORK, "candidates_genes.tsv"), sep="\t", index_col=0)
designs = pd.read_csv(os.path.join(WORK, "designs_raw.tsv"), sep="\t")
spec = pd.read_csv(os.path.join(RES, "specificity_report.tsv"), sep="\t")
pool = pd.read_csv(os.path.join(RES, "candidates_pool.tsv"), sep="\t")
final = pd.read_csv(os.path.join(RES, "final_multiplex.tsv"), sep="\t")
reps = pd.read_csv(os.path.join(WORK, "replicons.tsv"), sep="\t")
dnaa = pd.read_csv(os.path.join(WORK, "cds_table.tsv"), sep="\t")
dnaa = dnaa[dnaa.gene == "dnaA"].drop_duplicates("strain").set_index("strain")
S = sorted(final.strain)

# =============================================================================
# Figure 1 - workflow
# =============================================================================
either = (cds["class"] != "") | (cds.text_class != "")
eligible = either & (cds.replicon == "Chromosome") & ~cds.pseudo & ~cds.near_mobile & (cds.end - cds.start + 1 >= 300)
cand_primary = cand[cand.single_copy_nt & (cand.seqid == cand.strain.map(dnaa.seqid))]
uniq_pct = {}
for s in S:
    m = np.load(os.path.join(WORK, f"mask_{s}_{dnaa.seqid[s]}.npy"))
    uniq_pct[s] = 100 * m.mean()
n_combo = int(np.prod(pool.groupby("strain").size().values))

# worst heterodimer dG between oligos of two different assays, from the step-05 matrix
_dm = pd.read_csv(os.path.join(RES, "dimer_matrix.tsv"), sep="\t", index_col=0)
_other = np.array([[a.split("_")[0] != b.split("_")[0] for b in _dm.columns] for a in _dm.index])
worst_between_txt = f"{_dm.values[_other].min():.2f}".replace("-", "−")

steps = [
    ("7 complete genomes (PRJNA357031), RefSeq PGAP annotation",
     f"{len(cds):,} CDS  ·  {reps.role.eq('Chromosome').sum()} chromosomes, {reps.role.eq('Plasmid').sum()} plasmids", None),
    None,  # the two routes, drawn separately
    ("Eligible candidate genes",
     f"{eligible.sum():,} genes  ·  either route, chromosomal, not pseudo, ≥ 300 bp, > 5 kb from mobile elements", None),
    ("Sequence-aware filter: nucleotide uniqueness + single copy",
     f"{len(cand_primary):,} genes on the dnaA chromosome  ·  ≥ 150 bp unique DNA, single BLAST hit in own genome", None),
    ("Primer / probe design (primer3)",
     f"{len(designs):,} assays on {designs.locus_tag.nunique()} genes  ·  40 core + 40 unique genes per strain, 70–150 bp amplicons", None),
    ("Exhaustive in-silico PCR (bowtie1, ≤ 3 mismatches, 7 genomes + plasmids)",
     f"{int(spec['pass'].sum()):,} assays pass  ·  single on-target product, no near-perfect off-target primer site", None),
    ("Multiplex selection",
     f"{len(pool)} assays (top 6 genes per strain)  →  {n_combo:,} 7-way combinations scored on heterodimer ΔG", None),
    ("Final 7-plex, verified in silico",
     f"pooled in-silico PCR: exactly 7 products  ·  knock-out 7/7  ·  worst between-assay ΔG {worst_between_txt} kcal/mol", "final"),
]

fig, ax = plt.subplots(figsize=(7.2, 7.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 102); ax.axis("off")
W, H, X0 = 91, 7.6, 3
ys = [93, 76.5, 63, 51, 39, 27, 15, 3]


def box(x, y, w, h, title, sub, accent=None):
    ec = ORANGE if accent == "final" else MUTED
    lw = 1.6 if accent == "final" else 0.8
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.2",
                                fc=SURF, ec=ec, lw=lw))
    ax.text(x + 1.6, y + h * 0.66, title, fontsize=8.2, fontweight="bold", va="center", color=INK)
    ax.text(x + 1.6, y + h * 0.28, sub, fontsize=7, va="center", color=INK2, wrap=True)


def arrow(y_from, y_to, x=X0 + W / 2):
    ax.annotate("", xy=(x, y_to), xytext=(x, y_from),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9, shrinkA=0, shrinkB=0))


for i, st in enumerate(steps):
    if st is None:
        continue
    box(X0, ys[i], W, H, st[0], st[1], st[2])

# the two annotation routes + the k-mer track, side by side
rw, gap = (W - 2 * 2) / 3, 2
routes = [
    ("A · GFF text parsing\n(gene symbol / product)", [f"unique name/product: {int((cds.text_class == 'text_unique').sum()):,}",
                              f"named once in all 7: {int((cds.text_class == 'text_core_1x').sum()):,}"]),
    ("B · Protein homology\n(mmseqs2 all-vs-all)", [f"no homolog in other 6: {int((cds['class'] == 'A_unique').sum()):,}",
                                       f"1:1 in all 7 genomes: {int((cds['class'] == 'B_universal_1to1').sum()):,}"]),
    ("C · DNA uniqueness\n(18-mers, jellyfish)", ["unique fraction of primary",
                                           f"chromosome: {min(uniq_pct.values()):.0f}–{max(uniq_pct.values()):.0f} %"]),
]
ry, rh = ys[1] - 2.2, H + 4.4
for k, (t, lines) in enumerate(routes):
    x = X0 + k * (rw + gap)
    ax.add_patch(FancyBboxPatch((x, ry), rw, rh, boxstyle="round,pad=0.25,rounding_size=1.2",
                                fc="#f5f4f1", ec=GRID, lw=0.8))
    ax.text(x + 1.4, ry + rh - 1.2, t, fontsize=7.6, fontweight="bold", va="top", linespacing=1.15)
    for j, ln in enumerate(lines):
        ax.text(x + 1.4, ry + rh - 8.0 - 3.0 * j, ln, fontsize=7, color=INK2, va="center")
    ax.annotate("", xy=(x + rw / 2, ry + rh + 0.4), xytext=(X0 + W / 2, ys[0] - 0.3),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9, connectionstyle="arc3,rad=0"))
# A and B feed "eligible"; C feeds the sequence-aware filter
for k in (0, 1):
    x = X0 + k * (rw + gap) + rw / 2
    ax.annotate("", xy=(x, ys[2] + H + 0.3), xytext=(x, ry - 0.3),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9))
# C feeds the sequence-aware filter: dashed path around the right edge
xr, ym = X0 + W + 2.2, ys[3] + H / 2
ax.plot([X0 + W + 0.3, xr, xr], [ry + rh / 2, ry + rh / 2, ym], color=MUTED, lw=0.9, ls="--")
ax.annotate("", xy=(X0 + W + 0.4, ym), xytext=(xr, ym),
            arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9, ls="--"))
for i in range(2, len(ys) - 1):
    arrow(ys[i] - 0.3, ys[i + 1] + H + 0.3, x=X0 + W / 3)
save(fig, "fig1_workflow")

# =============================================================================
# Figure 2 - chromosome maps
# =============================================================================
BIN = 10_000
passing = spec[spec["pass"]]
fig = plt.figure(figsize=(7.2, 4.9))
for k, s in enumerate(S + ["legend"]):
    ax = fig.add_subplot(2, 4, k + 1, projection="polar")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_ylim(0, 1.0); ax.axis("off")
    if s == "legend":
        break
    seqid, ori = dnaa.seqid[s], dnaa.start[s]
    m = np.load(os.path.join(WORK, f"mask_{s}_{seqid}.npy"))
    L = len(m)
    theta = lambda pos: 2 * np.pi * (((pos - ori) % L) / L)
    # mid-replichore band (rel_ori 0.3-0.7) on both replichores
    for a0, a1 in ((0.15, 0.35), (0.65, 0.85)):
        ax.bar((a0 + a1) * np.pi, 0.28, width=(a1 - a0) * 2 * np.pi, bottom=0.54, color="#f2f1ee",
               edgecolor="none", zorder=0)
    # outer ring: fraction of unique positions per 10 kb
    nb = L // BIN
    frac = m[:nb * BIN].reshape(nb, BIN).mean(axis=1)
    centers = 2 * np.pi * (np.arange(nb) * BIN + BIN / 2) / L  # bins start at position 1
    centers = (centers - 2 * np.pi * (ori / L)) % (2 * np.pi)
    ax.bar(centers, 0.13, width=2 * np.pi * BIN / L, bottom=0.86, color=SEQ(frac), edgecolor="none",
           linewidth=0, zorder=2)
    # passing assays (ticks) on this chromosome
    pa = passing[(passing.strain == s) & (passing.seqid == seqid)]
    for pos in pa.amp_start:
        ax.plot([theta(pos)] * 2, [0.60, 0.80], color=GREY, lw=0.5, zorder=3)
    # chosen target: bold tick across the assay ring + dot on the outer ring
    f = final[final.strain == s].iloc[0]
    tt = theta(f.amp_start)
    ax.plot([tt, tt], [0.58, 0.84], color=ORANGE, lw=2.2, zorder=5, solid_capstyle="round")
    ax.scatter([tt], [0.84], s=28, color=ORANGE, zorder=6, edgecolors=SURF, linewidths=0.9)
    # ori / ter
    ax.scatter([0], [0.99], marker="v", s=20, color=INK, zorder=6)
    ax.text(0, 1.13, "ori", ha="center", va="center", fontsize=6.3, color=INK)
    ax.text(np.pi, 1.12, "ter", ha="center", va="center", fontsize=6.3, color=INK2)
    gene = f.gene if isinstance(f.gene, str) else f.locus_tag.split("_")[-1]
    sp = species[s].split()
    ax.text(0.5, 0.60, s, transform=ax.transAxes, ha="center", va="center", fontsize=8.5, fontweight="bold")
    ax.text(0.5, 0.50, f"{sp[0][0]}. {sp[1]}", transform=ax.transAxes, ha="center", va="center",
            fontsize=6.3, style="italic", color=INK2)
    ax.text(0.5, 0.40, f"{L / 1e6:.2f} Mb", transform=ax.transAxes, ha="center", va="center",
            fontsize=5.8, color=MUTED)
    ax.text(0.5, -0.16, f"target {gene} · rel_ori {f.rel_ori:.2f}", transform=ax.transAxes,
            ha="center", fontsize=6.6, color=INK)
# legend panel
lax = fig.axes[-1]
lax.remove()
lax = fig.add_axes([0.765, 0.08, 0.2, 0.36])
lax.axis("off")
cb = lax.inset_axes([0.05, 0.86, 0.85, 0.07])
cb.imshow(np.linspace(0, 1, 256)[None, :], cmap=SEQ, aspect="auto")
cb.set_yticks([]); cb.set_xticks([0, 255]); cb.set_xticklabels(["0", "1"], fontsize=6.5)
cb.tick_params(length=2)
lax.text(0.05, 1.0, "Outer ring: unique fraction\n(18-mers, 10 kb bins)", fontsize=6.8, va="bottom",
         transform=lax.transAxes)
lax.plot([0.05, 0.05], [0.55, 0.68], color=GREY, lw=1, transform=lax.transAxes)
lax.text(0.12, 0.615, "assay passing specificity", fontsize=6.8, va="center", transform=lax.transAxes)
lax.plot([0.05, 0.05], [0.38, 0.50], color=ORANGE, lw=1.6, transform=lax.transAxes)
lax.scatter([0.05], [0.50], s=22, color=ORANGE, transform=lax.transAxes, edgecolors=SURF, linewidths=0.8)
lax.text(0.12, 0.44, "chosen target\n(label: gene, rel_ori)", fontsize=6.8, va="center", transform=lax.transAxes)
lax.add_patch(Rectangle((0.015, 0.22), 0.07, 0.1, color="#f2f1ee", transform=lax.transAxes))
lax.text(0.12, 0.27, "mid-replichore\n(rel_ori 0.3–0.7)", fontsize=6.8, va="center", transform=lax.transAxes)
lax.scatter([0.05], [0.08], marker="v", s=22, color=INK, transform=lax.transAxes)
lax.text(0.12, 0.08, "dnaA (ori) at top", fontsize=6.8, va="center", transform=lax.transAxes)
fig.subplots_adjust(wspace=0.35, hspace=0.45, left=0.02, right=0.98, top=0.94, bottom=0.08)
save(fig, "fig2_chromosome_maps")

# =============================================================================
# Figure 3 - specificity
# =============================================================================
b = pd.read_csv(os.path.join(WORK, "final.bowtie.tsv"), sep="\t", header=None,
                names=["name", "strand", "chrom", "start", "seq", "qual", "n", "mism"], keep_default_na=False)
b["mm"] = b.mism.apply(lambda x: len(x.split(",")) if x else 0)
b["assay"] = b.name.str.split("_").str[0]
b["genome"] = b.chrom.str.split("|").str[0]
best = b.groupby(["assay", "genome"]).mm.min().unstack().reindex(index=S, columns=S)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw=dict(wspace=0.45))
for ax in (a1, a2):
    ax.set_xlim(-0.5, 6.5); ax.set_ylim(6.5, -0.5); ax.set_aspect("equal")
    ax.set_xticks(range(7)); ax.set_yticks(range(7))
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
amp_len = final.set_index("strain").amp_len
for i, a in enumerate(S):
    for j, g in enumerate(S):
        if i == j:
            a1.add_patch(Rectangle((j - 0.46, i - 0.46), 0.92, 0.92, fc=BLUE, ec="none"))
            a1.text(j, i, f"{amp_len[a]}\nbp", ha="center", va="center", fontsize=6.2, color="white", linespacing=0.95,
                    fontweight="bold")
        else:
            a1.add_patch(Rectangle((j - 0.46, i - 0.46), 0.92, 0.92, fc="#f5f4f1", ec="none"))
            v = best.loc[a, g]
            a1.text(j, i, "≥4" if pd.isna(v) else f"{int(v)}", ha="center", va="center", fontsize=6.8,
                    color=INK2)
a1.set_xticklabels(S); a1.set_yticklabels([f"{s} assay" for s in S])
a1.set_xlabel("genome (all replicons)")
a1.set_title("In-silico PCR per assay", loc="left", fontsize=8.5)
panel_letter(a1, "A", x=-0.22)

# B: knock-out - mask the target amplicon of the row genome, which assays still amplify?
sites = b.copy()
sites["end"] = sites.start + sites.seq.str.len()
ko = np.zeros((7, 7), dtype=bool)
fx = final.set_index("strain")
for i, s in enumerate(S):
    r = fx.loc[s]
    hit = (sites.chrom == f"{s}|{r.seqid}") & (sites.end > r.amp_start - 1) & (sites.start < r.amp_end)
    kept = sites[~hit]
    for j, a in enumerate(S):
        own = kept[kept.assay == a]
        f_ok = ((own.name == f"{a}_fwd") & (own.chrom == f"{a}|{fx.loc[a].seqid}") & (own.mm == 0)).any()
        r_ok = ((own.name == f"{a}_rev") & (own.chrom == f"{a}|{fx.loc[a].seqid}") & (own.mm == 0)).any()
        ko[i, j] = f_ok and r_ok
for i in range(7):
    for j in range(7):
        if ko[i, j]:
            a2.add_patch(Rectangle((j - 0.46, i - 0.46), 0.92, 0.92, fc=BLUE, ec="none"))
        else:
            a2.add_patch(Rectangle((j - 0.46, i - 0.46), 0.92, 0.92, fc=SURF, ec=ORANGE, lw=1.2, hatch="////"))
            a2.text(j, i, "lost", ha="center", va="center", fontsize=6.4, color=ORANGE, fontweight="bold",
                    bbox=dict(fc=SURF, ec="none", pad=0.6))
a2.set_xticklabels(S); a2.set_yticklabels([f"{s} masked" for s in S])
a2.set_xlabel("assay")
a2.set_title("Knock-out test (pooled primers)", loc="left", fontsize=8.5)
panel_letter(a2, "B", x=-0.22)
fig.text(0.07, -0.02, "A: blue = in-silico PCR product (designed length); grey cells give the fewest mismatches of any "
         "forward/reverse primer site in that genome (bowtie1 exhaustive search, ≥4 = no site with ≤3). "
         "No off-target product with ≤3 mismatches per primer.\nB: each row masks one target amplicon; "
         "blue = product still formed; only the masked target is lost.",
         fontsize=6.4, color=INK2, va="top", wrap=True)
save(fig, "fig3_specificity")

# =============================================================================
# Figure 4 - multiplex compatibility
# =============================================================================
order = list(final.sort_values(["well_4plus3", "strain"]).strain)
ol = [(f"{s} {role[0].upper() if role != 'probe' else 'P'}", fx.loc[s][role], s)
      for s in order for role in ("fwd", "rev", "probe")]
n = len(ol)
M = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        M[i, j] = primer3.calc_heterodimer(ol[i][1], ol[j][1], **COND).dg / 1000.0

# all combinations: worst between-assay dG (same definition as step 05)
oligos = {i: [r.fwd, r.rev, r.probe] for i, r in pool.iterrows()}
cache = {}


def dg(a, c):
    k = (a, c) if a < c else (c, a)
    if k not in cache:
        cache[k] = primer3.calc_heterodimer(a, c, **COND).dg / 1000.0
    return cache[k]


npool = len(pool)
ANY = np.zeros((npool, npool))
for i, j in itertools.combinations(range(npool), 2):
    if pool.strain[i] != pool.strain[j]:
        ANY[i, j] = ANY[j, i] = min(dg(x, y) for x in oligos[i] for y in oligos[j])
groups = [list(pool.index[pool.strain == s]) for s in S]
worst = np.array([ANY[np.ix_(c, c)].min() for c in map(list, itertools.product(*groups))])
chosen_idx = [pool.index[pool.set_id == sid][0] for sid in final.set_id]
chosen_val = ANY[np.ix_(chosen_idx, chosen_idx)].min()

fig = plt.figure(figsize=(7.2, 3.9))
a1 = fig.add_axes([0.07, 0.12, 0.47, 0.8])
a2 = fig.add_axes([0.75, 0.2, 0.24, 0.62])
vmin = -10
im = a1.imshow(np.clip(-M, 0, -vmin), cmap=SEQ, vmin=0, vmax=-vmin)
for k in range(1, 7):
    lw, c = (1.6, INK) if k == 4 else (0.8, SURF)
    a1.axhline(3 * k - 0.5, color=c, lw=lw); a1.axvline(3 * k - 0.5, color=c, lw=lw)
for k in range(7):
    a1.add_patch(Rectangle((3 * k - 0.5, 3 * k - 0.5), 3, 3, fill=False, ec=INK2, lw=0.9, ls=(0, (1.5, 1.5))))
a1.set_xticks(range(n)); a1.set_yticks(range(n))
a1.set_xticklabels([o[0] for o in ol], rotation=90, fontsize=5.8)
a1.set_yticklabels([o[0] for o in ol], fontsize=5.8)
a1.tick_params(length=0)
for sp in a1.spines.values():
    sp.set_visible(False)
a1.text(5.5, -1.4, "well W1", ha="center", fontsize=7, fontweight="bold")
a1.text(16.5, -1.4, "well W2", ha="center", fontsize=7, fontweight="bold")
# annotate the worst between-assay pair
between = np.array([[ol[i][2] != ol[j][2] for j in range(n)] for i in range(n)])
wi, wj = np.unravel_index(np.where(between, M, 0).argmin(), M.shape)
a1.add_patch(Rectangle((wj - 0.5, wi - 0.5), 1, 1, fill=False, ec=ORANGE, lw=1.4))
cax = fig.add_axes([0.555, 0.12, 0.012, 0.8])
cbar = fig.colorbar(im, cax=cax)
cbar.set_ticks([0, 2, 4, 6, 8, 10]); cbar.set_ticklabels(["0", "−2", "−4", "−6", "−8", "≤−10"])
cbar.set_label("heterodimer ΔG (kcal/mol, 37 °C)", fontsize=7)
cbar.outline.set_visible(False)
panel_letter(a1, "A", x=-0.02, y=1.045)
a1.text(0.0, -0.19, f"orange box: strongest between-assay pair ({ol[wi][0]} × {ol[wj][0]}, "
        + f"{M[wi, wj]:.2f}".replace("-", "−") + " kcal/mol); "
        + "dotted squares: within-assay oligos", transform=a1.transAxes, fontsize=6.2, color=INK2)

bins = np.arange(np.floor(worst.min()), 0.01, 0.25)
a2.hist(worst, bins=bins, color=GREY, edgecolor=SURF, linewidth=0.5)
a2.axvline(chosen_val, color=ORANGE, lw=1.6)
a2.text(chosen_val, a2.get_ylim()[1] * 0.97, f" chosen\n {chosen_val:.2f}".replace("-", "−"), color=ORANGE, fontsize=6.8,
        va="top", ha="left" if chosen_val < -3 else "right")
a2.set_xlabel("worst between-assay ΔG\nin the 7-plex (kcal/mol)")
a2.set_ylabel("combinations")
a2.spines[["top", "right"]].set_visible(False)
a2.grid(axis="y", color=GRID, lw=0.6); a2.set_axisbelow(True)
a2.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1000:.0f}k" if v else "0"))
# the add-on set (results folder with candidates_pool_wide.tsv) was chosen from a much wider
# pool than the one enumerated here, so the title must not claim "all" candidates
if os.path.exists(os.path.join(RES, "candidates_pool_wide.tsv")):
    a2.set_title(f"Reduced pool only: {len(worst):,}\n7-plexes (6 genes per strain)", loc="left", fontsize=8.5)
else:
    a2.set_title(f"All {len(worst):,}\ncandidate 7-plexes", loc="left", fontsize=8.5)
panel_letter(a2, "B", x=-0.2, y=1.1)
save(fig, "fig4_multiplex_compatibility")

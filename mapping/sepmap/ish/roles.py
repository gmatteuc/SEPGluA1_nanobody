"""Subunit genes against localisation genes: does the map need the second set?

Whether the nano map tracks GluA1 at the membrane or GluA1 as such:

    A    nano is proportional to total receptor
    B    nano is proportional to surface receptor = total x surface fraction

Both predict a postsynaptic map, so the postsynaptic-against-presynaptic contrast
of ish.words cannot separate them; a contrast within the postsynaptic compartment
can. The receptor's own subunit genes (Gria1-4) set its abundance; TARPs,
cornichons, CKAMP44, the PSD scaffolds that anchor the receptor and NSF set where
it goes. Under A the subunits carry the map and the localisation genes add nothing;
under B the localisation genes explain variance the subunits do not. It is a
statement about sets, not genes: one gene (Cacng8) can outrank another (Gria1) for
any number of reasons.

The roles (ROLES) are curated because the `category` column cannot express this:
it puts the metabotropic receptors Grm1-5 beside the TARPs and fills "trafficking"
with presynaptic vesicle machinery. They follow what each protein does and no rho
was consulted, though they were written after the per-gene ranking had been seen;
the permutation guards against that. Three steps, each worth more than the last:

    1. rho by role, the descriptive picture;
    2. commonality: the map rank-regressed on the subunit composite, on the
       localisation composite and on both, the explained variance split into what
       each explains alone and what they share; the composites are co-expressed,
       so the shared part is large by construction and the unique parts inform;
    3. a within-family permutation: the 19 family genes split 4 against 15 in all
       3,876 ways, which holds fixed the co-expression that makes the usual
       gene-category null useless (Fulcher 2021); exact, not sampled.

A sensitivity run (4) repeats 2 and 3 without the barely expressed family genes.
Presynaptic vesicle machinery is the specificity control: if the map liked
membrane trafficking of any kind, it would score too. Gria1 mRNA is not GluA1
protein, so a subunit set explaining less than localisation also fits a regional
translation gradient; only a real total-receptor channel separates the two, and
the green channel is not one (adult.sep_channel_check).

Run by run_ish_roles.py.
"""

import csv
import itertools
import os
from collections import defaultdict

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr

from sepmap.config import DATA

NANO = os.path.join(
    DATA, "comparisons_v2", "young_vs_adult", "region_means_per_mouse.csv"
)
GENES = os.path.join(DATA, "adult_v2", "ish", "gene_region_table.csv")
OUT = os.path.join(DATA, "adult_v2", "ish")

# the adult groups, pooled
ADULT_GROUPS = ("naive", "rws")

# the reading the project has settled on, for the tests and the figure
READING = "zref"
ALL_READINGS = ("zref", "cref", "subref", "ratio", "sepratio")

# 200 um voxels a structure's gene value needs
MIN_ISH_VOXELS = 10

# curated roles, from what the protein does, not from any correlation: subunit and
# localisation carry the argument, presyn is the specificity control, the rest are
# context and enter no test
ROLES = {
    # the receptor itself: these set how much receptor there is
    "Gria1": "subunit",
    "Gria2": "subunit",
    "Gria3": "subunit",
    "Gria4": "subunit",
    # auxiliary subunits: TARPs, cornichons, CKAMP44. They do not make receptor,
    # they decide how much of it reaches and stays in the synaptic membrane
    "Cacng3": "localisation",
    "Cacng4": "localisation",
    "Cacng5": "localisation",
    "Cacng7": "localisation",
    "Cacng8": "localisation",
    "Cnih2": "localisation",
    "Cnih3": "localisation",
    "Shisa9": "localisation",
    # the PSD scaffolds that anchor AMPARs, and NSF/synaptojanin, which hold
    # them at the surface and retrieve them
    "Dlg4": "localisation",
    "Dlg2": "localisation",
    "Grip1": "localisation",
    "Grip2": "localisation",
    "Homer1": "localisation",
    "Nsf": "localisation",
    "Synj1": "localisation",
    # vesicle fusion machinery: membrane trafficking, but presynaptic; the control
    # that separates AMPAR localisation from trafficking in general
    "Snap25": "presyn",
    "Stx1a": "presyn",
    "Syn1": "presyn",
    "Cplx1": "presyn",
    "Bsn": "presyn",
    "Vamp2": "presyn",
    "Syt1": "presyn",
    "Syt2": "presyn",
    "Syt3": "presyn",
    # trans-synaptic organisers
    "Nrxn1": "adhesion",
    "Lrfn2": "adhesion",
    "Ptprs": "adhesion",
    "Cbln2": "adhesion",
    "Nptx1": "adhesion",
    # other glutamate receptors: NMDA, kainate, metabotropic
    "Grin1": "other_glut_r",
    "Grin2a": "other_glut_r",
    "Grin2b": "other_glut_r",
    "Grin2c": "other_glut_r",
    "Grin2d": "other_glut_r",
    "Grin3a": "other_glut_r",
    "Grik1": "other_glut_r",
    "Grik2": "other_glut_r",
    "Grik3": "other_glut_r",
    "Grik4": "other_glut_r",
    "Grik5": "other_glut_r",
    "Grm1": "other_glut_r",
    "Grm2": "other_glut_r",
    "Grm3": "other_glut_r",
    "Grm4": "other_glut_r",
    "Grm5": "other_glut_r",
    # excitatory identity, immediate-early genes, and the original controls
    "Slc17a7": "excit_marker",
    "Slc17a6": "excit_marker",
    "Camk2a": "excit_marker",
    "Nrgn": "excit_marker",
    "Arc": "ieg",
    "Bdnf": "ieg",
    "Egr1": "ieg",
    "Fos": "ieg",
    "Nos1": "ieg",
}

# the roles of the 19 genes the permutation shuffles
FAMILY = ("subunit", "localisation")

# the roles in the order they are printed and drawn, and their labels
ORDER = [
    "subunit",
    "localisation",
    "presyn",
    "adhesion",
    "other_glut_r",
    "excit_marker",
    "ieg",
    "control_inhib",
    "control_glia",
    "control_struct",
]
NICE = {
    "subunit": "AMPAR subunits\nGria1-4",
    "localisation": "AMPAR localisation\nTARPs, PSD, NSF",
    "presyn": "presynaptic\nvesicle machinery",
    "adhesion": "trans-synaptic\norganisers",
    "other_glut_r": "other glutamate\nreceptors",
    "excit_marker": "excitatory\nidentity",
    "ieg": "immediate-early",
    "control_inhib": "inhibitory\ncontrols",
    "control_glia": "glial\ncontrols",
    "control_struct": "structural\ncontrols",
}


def adult_profile(reading):
    """{structure: mean over the adults} of one reading, from the per-mouse table.

    A structure's mean is over the adults that have it, however many they are.
    """
    per = defaultdict(list)
    with open(NANO, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["group"] in ADULT_GROUPS and r["reading"] == reading:
                per[r["structure"]].append(float(r["log2_value"]))
    return {s: float(np.mean(v)) for s, v in per.items()}


def gene_profiles():
    """{gene: {structure: expression}} and {gene: category}.

    Only structures covered by at least MIN_ISH_VOXELS voxels of the gene's grid.
    """
    out, cat = defaultdict(dict), {}
    with open(GENES, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if int(r["n_voxels"]) >= MIN_ISH_VOXELS:
                out[r["symbol"]][r["structure"]] = float(r["ish_mean"])
                cat[r["symbol"]] = r["category"]
    return out, cat


def role_of(gene, category):
    """Curated role where there is one, the original control label otherwise."""
    if gene in ROLES:
        return ROLES[gene]
    c = category.get(gene, "")
    return c if c.startswith("control") else "other"


def composite(genes, expr, structures):
    """One profile from a set of genes: the mean of their rank profiles."""
    r = [rankdata([expr[g][s] for s in structures]) for g in genes]
    return np.mean(r, axis=0)


def r2(y, xs):
    """Fraction of the variance of `y` explained linearly by the columns `xs`.

    Called on ranks, so it is a rank regression.
    """
    a = np.column_stack(list(xs) + [np.ones_like(y)])
    resid = y - a @ np.linalg.lstsq(a, y, rcond=None)[0]
    return float(1 - resid.var() / y.var())


def commonality(y, s, m):
    """R2 of each composite and of both, the part unique to each, the part shared.

    `s` is the subunit composite and `m` the localisation composite.
    """
    rs, rm, both = r2(y, [s]), r2(y, [m]), r2(y, [s, m])
    return dict(
        r2_subunit=rs,
        r2_localisation=rm,
        r2_both=both,
        unique_subunit=both - rm,
        unique_localisation=both - rs,
        shared=rs + rm - both,
    )


def permutation(y, expr, structures, family_genes, n_subunit, observed):
    """The statistic for every way of calling `n_subunit` family genes the subunit set.

    Splitting the same genes holds their co-expression fixed, the thing that makes
    the usual gene-category null meaningless here: each surrogate set is as
    co-expressed as the real one. The statistic is unique localisation minus
    unique subunit variance; p is the share of splits at or above `observed`, one
    added to both counts.
    """
    stats = []
    for pick in itertools.combinations(sorted(family_genes), n_subunit):
        rest = [g for g in family_genes if g not in pick]
        c = commonality(
            y, composite(pick, expr, structures), composite(rest, expr, structures)
        )
        stats.append(c["unique_localisation"] - c["unique_subunit"])
    stats = np.array(stats)
    p = float((np.sum(stats >= observed) + 1) / (len(stats) + 1))
    return stats, p


def sensitivity(y, expr, structures, family, roles):
    """Print the commonality and the permutation again, without the quiet family genes.

    The obvious objection to a negative result is that the localisation set carries
    passengers: Cacng5, Cacng7 and Grip2 are close to absent from the forebrain, so
    they can only dilute it. Dropping genes because they scored badly is the
    post-hoc move the permutation exists to catch, so the filter is on expression
    only: a gene's median energy across structures against the family's median,
    applied to both sets alike, never looking at a rho. If the conclusion survives
    this, it is not an artefact of dead weight.
    """
    # the family genes at or above the family's median expression
    level = {g: float(np.median([expr[g][t] for t in structures])) for g in family}
    cut = float(np.median(list(level.values())))
    kept = sorted(g for g in family if level[g] >= cut)
    sub = [g for g in kept if roles[g] == "subunit"]
    loc = [g for g in kept if roles[g] == "localisation"]
    print(
        f"\n4. sensitivity: family genes at or above the family median expression "
        f"({cut:.1f})"
    )
    print(
        f"   kept {len(kept)} of {len(family)}; dropped "
        f"{', '.join(g for g in family if g not in kept)}"
    )
    if len(sub) < 2 or len(loc) < 2:
        print("   too few genes left on one side -- not run")
        return
    print(f"   subunit ({len(sub)}) {', '.join(sub)}")
    print(f"   localisation ({len(loc)}) {', '.join(loc)}")

    # the two steps on the genes kept
    c = commonality(y, composite(sub, expr, structures), composite(loc, expr, structures))
    observed = c["unique_localisation"] - c["unique_subunit"]
    stats, p = permutation(y, expr, structures, kept, len(sub), observed)
    print(
        f"   unique subunit {c['unique_subunit']:+.3f}, "
        f"unique localisation {c['unique_localisation']:+.3f}, "
        f"shared {c['shared']:+.3f}"
    )
    print(f"   observed {observed:+.4f} against {len(stats):,} splits, p = {p:.4f}")


def main():
    """Run the three steps and the sensitivity run on the map, write, print and draw."""
    # the map, the gene profiles and each gene's role
    nano = adult_profile(READING)
    expr, category = gene_profiles()
    roles = {g: role_of(g, category) for g in expr}

    # every test uses the same structures: those the map and all 19 family
    # genes have in common, so nothing moves between comparisons
    family = sorted(g for g in expr if roles[g] in FAMILY)
    structures = sorted(set(nano).intersection(*[set(expr[g]) for g in family]))
    print(
        f"{len(expr)} genes, {len(family)} in the AMPAR family, "
        f"{len(structures)} shared structures, reading {READING}"
    )
    sub = sorted(g for g in family if roles[g] == "subunit")
    loc = sorted(g for g in family if roles[g] == "localisation")
    print(f"  subunit      ({len(sub):2d})  {', '.join(sub)}")
    print(f"  localisation ({len(loc):2d})  {', '.join(loc)}")

    # rho by role (step 1), every gene against every reading
    rows = []
    for reading in ALL_READINGS:
        prof = adult_profile(reading)
        for gene in sorted(expr):
            common = sorted(set(prof) & set(expr[gene]))
            if len(common) < 50:
                continue
            rho, _ = spearmanr([prof[s] for s in common], [expr[gene][s] for s in common])
            rows.append(
                dict(
                    reading=reading,
                    role=roles[gene],
                    symbol=gene,
                    n_structures=len(common),
                    rho=float(rho),
                )
            )
    path = os.path.join(OUT, "role_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(
            {k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()}
            for r in rows
        )
    print(f"{len(rows)} rows -> {path}")

    # the curated assignment, so it can be argued with
    with open(
        os.path.join(OUT, "gene_roles.csv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh)
        w.writerow(["symbol", "role", "original_category"])
        for g in sorted(expr):
            w.writerow([g, roles[g], category.get(g, "")])

    # step 1 printed for the reading of the tests
    here = [r for r in rows if r["reading"] == READING]
    by_role = defaultdict(list)
    for r in here:
        by_role[r["role"]].append(r["rho"])
    print(f"\n1. rho with the {READING} map, by curated role")
    for role in ORDER:
        v = by_role.get(role, [])
        if v:
            print(
                f"  {role:15s} n={len(v):3d}   median {np.median(v):+.3f}   "
                f"[{np.percentile(v, 25):+.3f} {np.percentile(v, 75):+.3f}]"
            )

    # the commonality of the two composites (step 2)
    y = rankdata([nano[s] for s in structures])
    s_comp = composite(sub, expr, structures)
    m_comp = composite(loc, expr, structures)
    c = commonality(y, s_comp, m_comp)
    print(
        f"\n2. variance of the {READING} map explained "
        f"(ranks, {len(structures)} structures)"
    )
    print(f"  subunit composite alone        R2 = {c['r2_subunit']:.3f}")
    print(f"  localisation composite alone   R2 = {c['r2_localisation']:.3f}")
    print(f"  both                           R2 = {c['r2_both']:.3f}")
    print(f"  unique to subunits             {c['unique_subunit']:+.3f}")
    print(f"  unique to localisation         {c['unique_localisation']:+.3f}")
    print(f"  shared                         {c['shared']:+.3f}")
    print(
        f"  composites correlate with each other at rho "
        f"{spearmanr(s_comp, m_comp).statistic:+.3f}"
    )

    # the within-family permutation (step 3)
    observed = c["unique_localisation"] - c["unique_subunit"]
    stats, p = permutation(y, expr, structures, family, len(sub), observed)
    print(
        f"\n3. every 4-of-{len(family)} split of the same family: {len(stats):,} of them"
    )
    print(f"  observed difference in unique variance {observed:+.4f}")
    print(
        f"  surrogate splits median {np.median(stats):+.4f}, "
        f"p95 {np.percentile(stats, 95):+.4f}"
    )
    print(f"  p = {p:.4f}  (fraction of splits reaching the functional one)")

    # the specificity control
    pre = by_role.get("presyn", [])
    print(
        f"\n  specificity control -- presynaptic vesicle machinery: "
        f"n={len(pre)}, median rho {np.median(pre):+.3f} "
        f"against localisation {np.median(by_role['localisation']):+.3f}"
    )

    # the sensitivity run (step 4) and the figure
    sensitivity(y, expr, structures, family, roles)
    figure(by_role, c, stats, observed, p, s_comp, m_comp, y)


def figure(by_role, c, stats, observed, p, s_comp, m_comp, y):
    """Draw steps 1 to 3 for the reading of the tests; saved as ish_roles.png."""
    fig, axes = plt.subplots(
        1, 3, figsize=(15.5, 5.2), gridspec_kw=dict(width_ratios=[1.7, 0.8, 1.0])
    )
    rng = np.random.default_rng(0)

    # left: rho by role, subunit genes in blue, localisation genes in red
    ax = axes[0]
    roles = [r for r in ORDER if by_role.get(r)]
    for i, role in enumerate(roles):
        v = by_role[role]
        colour = {"subunit": "#1f3b73", "localisation": "#c0392b"}.get(role, "0.65")
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.14, 0.14, len(v)),
            v,
            s=17,
            facecolor=colour,
            edgecolor="0.25",
            linewidth=0.4,
            zorder=2,
        )
        ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color="0.15", lw=1.8, zorder=3)
    ax.axhline(0, color="0.85", lw=0.7, zorder=0)
    ax.set_xticks(range(len(roles)))
    ax.set_xticklabels(
        [f"{NICE.get(r, r).replace(chr(10), ' ')} ({len(by_role[r])})" for r in roles],
        fontsize=7,
        rotation=32,
        ha="right",
    )
    ax.set_ylabel(f"Spearman with the adult {READING} map", fontsize=8)
    ax.set_title(
        "1. the map by curated role\n"
        "blue = what sets abundance, red = what sets localisation",
        fontsize=9,
    )

    # middle: the commonality
    ax = axes[1]
    parts = [c["unique_subunit"], c["shared"], c["unique_localisation"]]
    ax.bar(
        range(3),
        parts,
        color=["#1f3b73", "0.75", "#c0392b"],
        edgecolor="0.25",
        linewidth=0.5,
    )
    ax.set_xticks(range(3))
    ax.set_xticklabels(
        ["unique to\nsubunits", "shared", "unique to\nlocalisation"], fontsize=7.5
    )
    ax.set_ylabel("variance of the map explained (R2, ranks)", fontsize=8)
    for i, v in enumerate(parts):
        ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_title(
        f"2. commonality\ncomposites correlate at rho "
        f"{spearmanr(s_comp, m_comp).statistic:+.2f}",
        fontsize=9,
    )

    # right: the splits of the family, and the split by function
    ax = axes[2]
    ax.hist(stats, bins=60, color="0.72", edgecolor="0.35", linewidth=0.3)
    ax.axvline(observed, color="#c0392b", lw=2)
    ax.annotate(
        "the split\nby function",
        (observed, ax.get_ylim()[1] * 0.9),
        color="#c0392b",
        fontsize=7.5,
        ha="right",
        va="top",
        xytext=(-6, 0),
        textcoords="offset points",
    )
    ax.set_xlabel("unique(localisation) - unique(subunit)", fontsize=8)
    ax.set_ylabel(f"splits of the same {len(stats):,}", fontsize=8)
    ax.set_title(f"3. every 4-of-19 split of the family\np = {p:.4f}", fontsize=9)

    for ax in axes:
        ax.tick_params(labelsize=7)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    # the title says what the panels show
    fig.suptitle(
        "Does the map need the localisation genes, or do the subunits account for it?"
        "\nThe direction is as predicted -- the subunits add almost nothing beyond "
        "localisation -- but splitting the same 19 genes any other way does as well "
        "(panel 3), so this panel is no evidence either way.",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    path = os.path.join(OUT, "ish_roles.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"\n{path}")

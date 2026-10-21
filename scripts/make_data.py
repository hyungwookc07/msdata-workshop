"""Hands-on data from Small et al., J Lipid Res 2024 (PXD054323; github.com/lewinsmall/Lipogenesis-Proteomics).

4-week cohort, liver, MaxQuant LFQ.
Steps: remove reverse / contaminants / only-by-site, unique peptides >= 2, LFQ 0 -> NA, log2,
median normalization, flag the 5 muscle-contaminated samples (authors' list), keep proteins with
>= 70% valid values in both compared groups, limma-style moderated t (two-group, NA-aware, no imputation),
BH-FDR. Comparison: high-fat (HFD) vs high-starch (HSD).
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats, special

PG, META, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(f"{OUT}/data", exist_ok=True)

pg = pd.read_csv(PG, sep="\t", low_memory=False)
pg = pg[(pg["Reverse"] != "+") & (pg["Potential contaminant"] != "+") & (pg["Only identified by site"] != "+")]
pg = pg[pg["Unique peptides"] >= 2].copy()
pg["gene"] = pg["Gene names"].fillna("").astype(str).str.split(";").str[0].replace("", np.nan)
pg["protein"] = pg["Majority protein IDs"].astype(str).str.split(";").str[0]

meta = pd.read_csv(META, encoding="utf-8-sig")
meta.columns = [c.strip() for c in meta.columns]
lfq_cols = {c: int(c.split("_")[-1]) for c in pg.columns if c.startswith("LFQ intensity")}
meta = meta[(meta.time_on_diet == "four_weeks") & meta.Animal_ID.isin(lfq_cols.values())].copy()
diet_label = {"fat": "HFD", "starch": "HSD", "chow": "Chow"}
meta["group"] = meta.Diet.map(diet_label)
outliers = [20842, 21167, 20833, 20821, 21182]          # muscle contamination (authors' script)
meta["muscle_contaminated"] = meta.Animal_ID.isin(outliers)
meta = meta.sort_values(["group", "Animal_ID"])
meta["sample"] = meta.group + "_" + meta.Animal_ID.astype(str)

inv = {v: k for k, v in lfq_cols.items()}
raw = pg[[inv[a] for a in meta.Animal_ID]].replace(0, np.nan).astype(float)
raw.columns = meta["sample"].values
raw.index = pg["protein"].values
L = np.log2(raw)
genes = pd.Series(pg["gene"].values, index=pg["protein"].values)

# median normalization (diff.median)
N = L - L.median() + L.median().median()

# save matrices (4-week, all three diets) and sample sheet
L.round(4).to_csv(f"{OUT}/data/protein_matrix_raw_log2.csv")
N.round(4).to_csv(f"{OUT}/data/protein_matrix_norm_log2.csv")
pd.DataFrame({"protein": genes.index, "gene": genes.values}).to_csv(f"{OUT}/data/protein_info.csv", index=False)
meta[["sample", "Animal_ID", "group", "Diet", "liver_triglyceride_content", "liver_lipogenesis",
      "muscle_contaminated"]].to_csv(f"{OUT}/data/samples.csv", index=False)

# ---- DE (limma-style), excluding contaminated samples ---------------------
keep = meta[~meta.muscle_contaminated]


def trigamma_inverse(x):
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if -dif / y < 1e-8:
            break
    return y


def fit_f_dist(x, df1):   # limma::fitFDist
    keepv = (x > 0) & np.isfinite(x) & (df1 > 0)
    x, df1 = x[keepv], df1[keepv]
    e = np.log(x) - special.digamma(df1 / 2) + np.log(df1 / 2)
    emean = e.mean()
    evar = ((e - emean) ** 2).sum() / (len(e) - 1) - special.polygamma(1, df1 / 2).mean()
    if evar > 0:
        d0 = 2 * trigamma_inverse(evar)
        return float(np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))), float(d0)
    return float(np.exp(emean)), np.inf


def limma_two_group(ga, gb):
    A = N[keep[keep.group == ga]["sample"]]
    B = N[keep[keep.group == gb]["sample"]]
    ok = (A.notna().mean(axis=1) >= 0.7) & (B.notna().mean(axis=1) >= 0.7)
    A, B = A[ok], B[ok]
    n1, n2 = A.notna().sum(axis=1), B.notna().sum(axis=1)
    m1, m2 = A.mean(axis=1), B.mean(axis=1)
    coef = m1 - m2
    ss = ((A.sub(m1, axis=0)) ** 2).sum(axis=1) + ((B.sub(m2, axis=0)) ** 2).sum(axis=1)
    df = (n1 + n2 - 2).astype(float)
    s2 = ss / df
    s20, d0 = fit_f_dist(s2.values, df.values)
    s2_post = (df * s2 + d0 * s20) / (df + d0)
    t = coef / (np.sqrt(s2_post) * np.sqrt(1 / n1 + 1 / n2))
    p = pd.Series(2 * stats.t.sf(np.abs(t), df + d0), index=A.index)
    o = np.argsort(p.values)
    r = p.values[o] * len(p) / (np.arange(len(p)) + 1)
    bh = np.empty(len(p)); bh[o] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    res = pd.DataFrame({"protein": A.index, "gene": genes.loc[A.index].values, "log2FC": coef.values,
                        "AveExpr": pd.concat([A, B], axis=1).mean(axis=1).values, "t": t.values,
                        "P.Value": p.values, "adj.P.Val": bh, f"n_{ga}": n1.values, f"n_{gb}": n2.values})
    print(f"{ga} vs {gb}: n={A.shape[1]} vs {B.shape[1]}, tested {len(res)}, prior s0^2={s20:.4f} d0={d0:.2f}")
    return res.sort_values("P.Value").reset_index(drop=True)


de = limma_two_group("HFD", "HSD")
de.round({"log2FC": 4, "AveExpr": 3, "t": 4}).to_csv(f"{OUT}/data/de_results_limma.csv", index=False)
for ga in ["HFD", "HSD"]:
    limma_two_group(ga, "Chow").round({"log2FC": 4, "AveExpr": 3, "t": 4}).to_csv(f"{OUT}/data/de_{ga}_vs_Chow.csv", index=False)
FDR_CUT, LFC_CUT = 0.05, np.log2(1.5)
up = de[(de["adj.P.Val"] < FDR_CUT) & (de.log2FC >= LFC_CUT)].gene.dropna().unique()
dn = de[(de["adj.P.Val"] < FDR_CUT) & (de.log2FC <= -LFC_CUT)].gene.dropna().unique()
bg = de.gene.dropna().unique()
for name, arr in [("up_genes", up), ("down_genes", dn), ("background", bg)]:
    open(f"{OUT}/data/{name}.txt", "w").write("\n".join(arr) + "\n")
rk = de.dropna(subset=["gene"]).assign(a=lambda d: d.t.abs()).sort_values("a", ascending=False)
rk = rk.drop_duplicates("gene").sort_values("t", ascending=False)
rk[["gene", "t"]].round(4).to_csv(f"{OUT}/data/ranked_list.rnk", sep="\t", header=False, index=False)

print("samples (4wk):", meta.group.value_counts().to_dict(), "| excluded:", len(outliers))
print("proteins after contaminant/peptide filter:", len(L), "| missing %:", round(float(L.isna().mean().mean() * 100), 1))
print(f"FDR<0.05 & |FC|>=1.5: up {len(up)}, down {len(dn)} | FDR<0.05 any: {(de['adj.P.Val']<0.05).sum()}")
print(de.head(12)[["gene", "log2FC", "t", "adj.P.Val"]].to_string())

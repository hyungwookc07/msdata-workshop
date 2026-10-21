"""GSEA용 마우스 GO Biological Process gene set 파일 만들기 (data/genesets/GO_BP_mouse.gmt).

출처
- gene set: MSigDB C5 GO:BP v7.5.1 (R 패키지 msigdbr 7.5.1에 포함된 데이터, github.com/cran/msigdbr)
- 사람 -> 마우스 ortholog: R 패키지 babelgene (github.com/cran/babelgene), 3개 이상 DB가 지지하는 ortholog
R의 msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")와 같은 방식입니다.

사용법
  git clone --depth 1 --branch 7.5.1 https://github.com/cran/msigdbr
  git clone --depth 1 https://github.com/cran/babelgene
  pip install pyreadr
  python make_go_bp_gmt.py msigdbr/R/sysdata.rda babelgene/R/sysdata.rda data/genesets/GO_BP_mouse.gmt
"""
MS, BG, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
m = pyreadr.read_r(MS); b = pyreadr.read_r(BG)
gs = m["msigdbr_genesets"]; gg = m["msigdbr_geneset_genes"]; genes = m["msigdbr_genes"]
gs = gs[(gs.gs_cat == "C5") & (gs.gs_subcat == "GO:BP")]
o = b["orthologs_df"]; o = o[(o.taxon_id == 10090) & o.symbol.notna()]
o = o[o.support_n >= 3]            # msigdbr default style: keep orthologs supported by >=3 sources
x = gg[gg.gs_id.isin(gs.gs_id)].merge(genes[["gene_id", "human_gene_symbol"]], on="gene_id")
x = x.merge(o[["human_symbol", "symbol"]], left_on="human_gene_symbol", right_on="human_symbol")
x = x.merge(gs[["gs_id", "gs_name", "gs_exact_source"]], on="gs_id")
def pretty(n):
    s = n.replace("GOBP_", "").replace("_", " ").lower()
    return s
with open(OUT, "w") as f:
    for (name, go), d in x.groupby(["gs_name", "gs_exact_source"]):
        g = sorted(set(d.symbol))
        if len(g) >= 5:
            f.write(f"{pretty(name)} ({go})\tMSigDB_C5_GOBP_v7.5.1\t" + "\t".join(g) + "\n")
print("terms", x.gs_id.nunique())

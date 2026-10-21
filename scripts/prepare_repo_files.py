"""저장소에 올릴 백업 결과(ORA, STRING)를 만드는 스크립트 (인터넷이 되는 Colab에서 1회 실행).
GSEA용 GO-BP gene set(data/genesets/GO_BP_mouse.gmt)은 scripts/make_go_bp_gmt.py로 이미 만들어 두었습니다.

실행 전: data/ 폴더에 de_results_limma.csv, up_genes.txt, down_genes.txt, background.txt 가 있어야 합니다.
결과: data/backup/ora_up.csv, ora_down.csv, string_down.tsv
"""
import os, io, requests
import pandas as pd
from gprofiler import GProfiler

os.makedirs("data/backup", exist_ok=True)

# 1) ORA 백업 (g:Profiler, background = 정량된 단백질 전체)
read = lambda p: [l.strip() for l in open(p) if l.strip()]
bg = read("data/background.txt")
gpr = GProfiler(return_dataframe=True, user_agent="bertis_workshop")
for label in ["up", "down"]:
    res = gpr.profile(organism="mmusculus", query=read(f"data/{label}_genes.txt"), background=bg,
                      domain_scope="custom", sources=["GO:BP", "KEGG", "REAC"], user_threshold=0.05)
    res.to_csv(f"data/backup/ora_{label}.csv", index=False)
    print(f"ORA {label}:", len(res), "terms")

# 2) STRING 백업 (감소 단백질, score >= 0.7)
r = requests.post("https://string-db.org/api/tsv/network",
                  data={"identifiers": "\r".join(read("data/down_genes.txt")), "species": 10090,
                        "required_score": 700, "caller_identity": "bertis_workshop"}, timeout=60)
r.raise_for_status()
open("data/backup/string_down.tsv", "w").write(r.text)
print("STRING edges:", len(pd.read_csv(io.StringIO(r.text), sep="\t")))

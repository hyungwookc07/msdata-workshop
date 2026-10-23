# Bioinformatics Analysis for Proteomics: Hands-on

Hands-on materials for *Tutorial: Bioinformatics analysis for proteomics*, part of the session *From sample preparation to data processing for proteomics* (AI Informatics Team, Bertis).
Starting from differential expression (DE) results of a mouse liver proteome, you will run ORA, GSEA and a STRING network analysis and compare the results with the published paper.

**Dataset**: high-fat diet (HFD) vs high-starch diet (HSD), mouse liver, label-free DDA (MaxLFQ). Small et al., *J Lipid Res* 2024.

**Slides**: [`slides/Bertis_Proteomics_Tutorial.pptx`](slides/Bertis_Proteomics_Tutorial.pptx)

## Quick start

Pick one track. Both use the same data and give the same results.

### Web track (browser only)

Download the files from the `data/` folder (click a file, then the download button) and use them in each web tool.

| Step | Web tool | Settings | File |
|---|---|---|---|
| 1. DE results | Excel / Google Sheets | Sort by `adj.P.Val` | `data/de_results_limma.csv`, `figures/volcano_HFD_vs_HSD.png` |
| 2. ORA | [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) | Organism: *Mus musculus*; Advanced options → Statistical domain scope: Custom | Query: `data/down_genes.txt` (or `up_genes.txt`); Background: `data/background.txt` |
| 3. GSEA | [WebGestalt](https://www.webgestalt.org) | Method: GSEA; Functional database: geneontology → Biological Process; gene set size 100–500 | `data/ranked_list.rnk` |
| 4. Network | [STRING](https://string-db.org) | Multiple proteins; *Mus musculus*; minimum interaction score: high (0.700) | `data/down_genes.txt` |

### Code track (Google Colab)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/hyungwookc07/msdata-workshop/blob/main/hands_on_diet_liver.ipynb)

1. Click the badge above. A Google account is required; nothing to install.
2. Run the cells from top to bottom with `Shift + Enter`. No coding needed.
3. To keep your results: `File > Save a copy in Drive`.

| Section | Content |
|---|---|
| 0. Setup | Install packages, load data |
| 1. DE results | Result table, cutoffs, volcano plot |
| 2. ORA | g:Profiler, with all quantified proteins as background |
| 3. GSEA | gseapy prerank on GO Biological Process, ranked by limma t |
| 4. Network | STRING (mouse, score ≥ 0.7), hub proteins |
| 5. Link to the paper | HFD vs chow and HSD vs chow, side by side |
| Appendix | Raw data: missing values, normalization, PCA |

## What you should see

| Analysis | HSD↑ (lower in HFD) | HFD↑ (higher in HFD) |
|---|---|---|
| DE (FDR < 0.05, ≥ 1.5-fold) | 23 proteins, incl. Acaca, Acly, Fasn, Scd1, Me1, G6pdx | 12 proteins, incl. Ech1, Eci1, Cyp4a10 |
| ORA | Fatty acid biosynthesis, pyruvate metabolism | Weak (short list) |
| GSEA (GO-BP) | Lipid biosynthetic process, carbohydrate metabolic process | Organic acid, small molecule and lipid catabolic process |
| STRING | Lipogenic enzymes form one module | |

These match Fig. 4A of the paper. At FDR < 0.05 without a fold-change cutoff, we find 204 proteins (paper: 200).

## Data

Small L, et al. Liver proteomics identifies a disconnect between proteins associated with de novo lipogenesis and triglyceride storage. *J Lipid Res*. 2024;65:100687. PRIDE: PXD054323.
MaxQuant output and the authors' R analysis: [lewinsmall/Lipogenesis-Proteomics](https://github.com/lewinsmall/Lipogenesis-Proteomics)

- Male C57BL/6J mice, 4 weeks on diet, liver, label-free DDA (MaxLFQ), MaxQuant
- Samples: chow, HSD and HFD at 4 weeks, 15 each. The 5 samples with muscle contamination (flagged by the authors) are excluded from DE.
- Comparison: HFD (n = 14) vs HSD (n = 13); log2FC = HFD / HSD

### Processing (`scripts/make_data.py`)

1. Remove reverse, potential contaminant and only-identified-by-site entries; keep ≥ 2 unique peptides
2. LFQ intensity 0 → missing; log2 transform
3. Median normalization
4. Keep proteins measured in ≥ 70% of samples in both groups
5. limma-style moderated t-test (empirical Bayes, no imputation), BH-FDR

The paper used a three-group limma model in R with imputation, so individual numbers differ slightly. The conclusions are the same.

### Files

| File | Content |
|---|---|
| `data/protein_matrix_raw_log2.csv` | Protein × sample log2 LFQ, before normalization (45 samples) |
| `data/protein_matrix_norm_log2.csv` | After median normalization |
| `data/samples.csv` | Sample sheet: diet, liver TG, lipogenesis rate, muscle contamination |
| `data/protein_info.csv` | Protein ID ↔ gene symbol |
| `data/de_results_limma.csv` | DE results, HFD vs HSD (log2FC, t, P.Value, adj.P.Val) |
| `data/de_HFD_vs_Chow.csv`, `data/de_HSD_vs_Chow.csv` | DE results vs chow (section 5) |
| `data/down_genes.txt`, `data/up_genes.txt` | Lower / higher in HFD (FDR < 0.05, ≥ 1.5-fold) |
| `data/background.txt` | ORA background: all tested proteins |
| `data/ranked_list.rnk` | GSEA ranked list: gene and limma t (two columns, no header) |
| `data/genesets/GO_BP_mouse.gmt` | GO Biological Process gene sets for mouse (MSigDB C5 GO:BP v7.5.1, mouse orthologs; built with `scripts/make_go_bp_gmt.py`) |
| `figures/` | Volcano plot (HFD vs HSD) and the vs-chow comparison |
| `slides/Bertis_Proteomics_Tutorial.pptx` | Tutorial slides: bioinformatics analysis for proteomics |

## Notes

- The notebook reads data from this repository, so it runs the same way for everyone.
- The ORA and STRING sections call the g:Profiler and STRING servers. If a server is busy, wait a moment and run the cell again, or use the web track.
- Plot labels in the notebook are in English because Colab has no Korean font by default.

## References

- Small L, et al. *J Lipid Res*. 2024;65:100687.
- Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43:e47.
- Subramanian A, et al. Gene set enrichment analysis. *PNAS*. 2005;102:15545–15550.
- Kolberg L, et al. g:Profiler: interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update). *Nucleic Acids Res*. 2023;51:W207–W212.
- Szklarczyk D, et al. The STRING database in 2023. *Nucleic Acids Res*. 2023;51:D638–D646.
- Fang Z, et al. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*. 2023;39:btac757.

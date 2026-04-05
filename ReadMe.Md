GSDS for functional enrichment analyses
================
Dimitri Meistermann
2026-04-05

GSDS (Gene Set Differential Scoring) is focused on functional
enrichment. The package has 3 purposes: - provide functions to manage
gene set databases - implement well known methods to perform functional
enrichment: over-representation analysis (ORA) and gene set enrichment
analysis (GSEA) - provide a new method to perform functional enrichment
based on differential scoring (GSDS)

The GSDS method computes an activation score for each gene set in each
sample using a PCA (similar to Pathway Level analysis of Gene Expression
or PLAGE). This gene score is used to identify gene sets that are
differentially activated between two conditions.

## Installation

In a R console:

    install.packages("devtools")
    devtools::install_github("https://github.com/DimitriMeistermann/GSDS", build_vignettes = FALSE)

For a manual installation of dependencies:

    install.packages("BiocManager")
    BiocManager::install( c("SummarizedExperiment","grid",
    "AnnotationDbi","biomaRt","BiocParallel","ComplexHeatmap",
    "circlize","coop","data.table","fgsea","ggplot2","GO.db",
    "irlba","KEGGREST",""matrixStats","RSQLite","ggrepel","labeling",
    "scales","sparseMatrixStats","RSQLite","S4Vectors","knitr",
    "rmarkdown","reactome.db","qualpalr")

The package is ready to load.

``` r
library(GSDS)
```

## Over representation analysis (ORA)

For any functional enrichment method, we have to first download several
gene set databases by using `getDBterms`. We will use here *kegg* and
the *Gene Ontology* (biological process). Note that you can use your own
database of gene sets.

``` r
#species specific package will be load/asked for installation
geneSetDB<-getDBsets(species = "Human", database = c("kegg","goBP"))
#> 
#> Downloading Human gene sets from kegg, goBP returning SYMBOL IDs
#> 'select()' returned 1:1 mapping between keys and columns
geneSetDB$kegg[1:2]
#> $`hsa00010 Glycolysis / Gluconeogenesis`
#>  [1] "AKR1A1"  "ADH1A"   "ADH1B"   "ADH1C"   "ADH4"    "ADH5"    "ADH6"    "GALM"    "ADH7"    "LDHAL6A" "DLAT"    "DLD"     "ENO1"    "ENO2"   
#> [15] "ENO3"    "ALDH2"   "ALDH3A1" "ALDH1B1" "FBP1"    "ALDH3B1" "ALDH3B2" "ALDH9A1" "ALDH3A2" "ALDOA"   "ALDOB"   "ALDOC"   "G6PC1"   "GAPDH"  
#> [29] "GAPDHS"  "GCK"     "GPI"     "HK1"     "HK2"     "HK3"     "ENO4"    "LDHA"    "LDHB"    "LDHC"    "PGAM4"   "ALDH7A1" "PCK1"    "PCK2"   
#> [43] "PDHA1"   "PDHA2"   "PDHB"    "PFKL"    "PFKM"    "PFKP"    "PGAM1"   "PGAM2"   "PGK1"    "PGK2"    "PGM1"    "PKLR"    "PKM"     "PGM2"   
#> [57] "ACSS2"   "G6PC2"   "BPGM"    "TPI1"    "HKDC1"   "ADPGK"   "ACSS1"   "FBP2"    "LDHAL6B" "G6PC3"   "MINPP1" 
#> 
#> $`hsa00020 Citrate cycle (TCA cycle)`
#>  [1] "CS"     "DLAT"   "DLD"    "DLST"   "FH"     "IDH1"   "IDH2"   "IDH3A"  "IDH3B"  "IDH3G"  "MDH1"   "MDH2"   "ACLY"   "ACO1"   "OGDH"   "ACO2"  
#> [17] "PC"     "PCK1"   "PCK2"   "PDHA1"  "PDHA2"  "PDHB"   "OGDHL"  "SDHA"   "SDHB"   "SDHC"   "SDHD"   "SUCLG2" "SUCLG1" "SUCLA2"
```

Then we can make a simple over representation analysis.

``` r
data("DEgenesPrime_Naive")
vectorIsDE<-DEgenesPrime_Naive$isDE!="NONE";names(vectorIsDE)<-rownames(DEgenesPrime_Naive)
resEnrich<-enrich.ora(vectorIsDE,DBsets = geneSetDB)
head(resEnrich[order(resEnrich$pval,decreasing = FALSE),])
#>                                                                                             term database size_db size_universe obsOverlap expectOverlap
#> goBP.GO:0045926 negative regulation of growth           GO:0045926 negative regulation of growth     goBP      18            16         11      2.768088
#> goBP.GO:0071276 cellular response to cadmium ion     GO:0071276 cellular response to cadmium ion     goBP      23            19         12      3.287104
#> goBP.GO:0010273 detoxification of copper ion             GO:0010273 detoxification of copper ion     goBP      16            12          9      2.076066
#> goBP.GO:0006882 intracellular zinc ion homeostasis GO:0006882 intracellular zinc ion homeostasis     goBP      35            31         15      5.363170
#> goBP.GO:0035904 aorta development                                   GO:0035904 aorta development     goBP      33            31         15      5.363170
#> goBP.GO:0005975 carbohydrate metabolic process         GO:0005975 carbohydrate metabolic process     goBP     176           141         42     24.393773
#>                                                    OEdeviation   log2OR         pval       padj
#> goBP.GO:0045926 negative regulation of growth       0.06321216 3.399468 7.578169e-06 0.04724374
#> goBP.GO:0071276 cellular response to cadmium ion    0.06690559 3.039860 1.061419e-05 0.04724374
#> goBP.GO:0010273 detoxification of copper ion        0.05316830 3.846146 1.824225e-05 0.05413082
#> goBP.GO:0006882 intracellular zinc ion homeostasis  0.07400040 2.169698 6.602616e-05 0.11755298
#> goBP.GO:0035904 aorta development                   0.07400040 2.169698 6.602616e-05 0.11755298
#> goBP.GO:0005975 carbohydrate metabolic process      0.13519673 1.030602 1.732695e-04 0.23589224
```

# functional class scoring (FCS)

The Functional class scoring family of functional enrichment is here
implemented by a wrapper of `fGSEA`. For any FCS method, we have to
build first a score at the gene level. I propose here a score for
Differentially Expressed genes based on inverse of p-value with the sign
of the Log2(Fold-change): \$(1-p) \$

``` r
fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
resEnrich<-enrich.fcs(fcsScore,DBsets = geneSetDB,returnGenes = TRUE) #return genes of the gene set in the dataframe
head(resEnrich[order(resEnrich$pval,decreasing = FALSE),])
#>                                                                                     term database size         pval         padj   log2err        ES
#> kegg.hsa00190 Oxidative phosphorylation               hsa00190 Oxidative phosphorylation     kegg  138 5.629917e-12 5.011753e-08 0.8870750 0.5095867
#> goBP.GO:0006412 translation                                       GO:0006412 translation     goBP  311 1.210738e-10 5.388995e-07 0.8266573 0.3497149
#> goBP.GO:0042254 ribosome biogenesis                       GO:0042254 ribosome biogenesis     goBP  122 3.870420e-10 1.148483e-06 0.8140358 0.4661679
#> goBP.GO:0006364 rRNA processing                               GO:0006364 rRNA processing     goBP  162 1.191541e-08 2.651774e-05 0.7477397 0.3886989
#> kegg.hsa05016 Huntington disease                             hsa05016 Huntington disease     kegg  308 1.506403e-08 2.682000e-05 0.7337620 0.3331159
#> goBP.GO:1902600 proton transmembrane transport GO:1902600 proton transmembrane transport     goBP  181 3.064544e-08 4.212344e-05 0.7195128 0.4108231
#>                                                     NES        genes
#> kegg.hsa00190 Oxidative phosphorylation        2.680077 NDUFC2-K....
#> goBP.GO:0006412 translation                    2.165859 AARS1, A....
#> goBP.GO:0042254 ribosome biogenesis            2.511697 BYSL, C1....
#> goBP.GO:0006364 rRNA processing                2.244804 BYSL, CD....
#> kegg.hsa05016 Huntington disease               2.052540 NDUFC2-K....
#> goBP.GO:1902600 proton transmembrane transport 2.264197 SLC25A4,....

# The result can be visualized as a volcano plot
volcanoPlot(resEnrich,effectSizeCol = "NES", adjPvalCol = "padj",labelCol = resEnrich$term)
```

<img src="ReadMe_files/figure-gfm/FCS-1.png" alt="" style="display: block; margin: auto;" />

We can export enrichment with the genes of each gene set in a tsv file.

    exportEnrich(resEnrich,"test/resEnrich.tsv")

## Gene Set Differential Activation (GSDS)

We can also perform an enrichment based on gene set differential
activation with *GSDS*.

GSDS was developed independently from
[GSVA](https://bioconductor.org/packages/release/bioc/html/GSVA.html)
but is has the same idea of working on a matrix of gene set scores. As
*GSVA*, the first step is to translate the expression matrix in an
activation score matrix where each column is a gene set and each row is
a sample.

The particularity of *GSDS* is to perform this step by taking the first
Principal Componant of the matrix $genesOfGeneSet \times samples$ from a
Principal Component Analysis (PCA). One PCA is done per gene set and the
obtained values are the activation score for this particular gene set.
Hence, one PCA is performed for each gene set.

The activation score matrix is then used to perform differential
analysis with a linear model between two conditions.

`computeActivationScore` and `GSDS` return a `SummarizedExperiment`
object, wich contain the activation score as an assay, and the metadata
of the pathway, including the contributions in `rowData(object)`.

``` r
data("bulkLogCounts")
data("sampleAnnot")
geneSetActivScore<-computeActivationScore(bulkLogCounts,DBsets = geneSetDB)
resGSDS<-GSDS(geneSetActivScore = geneSetActivScore,colData = sampleAnnot,contrast = c("culture_media","T2iLGO","KSR+FGF2"),DBsets = geneSetDB)
resGSDS.sel <- resGSDS[abs(rowData(resGSDS)$log2FoldChange) > 8 & rowData(resGSDS)$padj < 0.001,]
nrow(resGSDS.sel)
#> [1] 83
head(rowData(resGSDS.sel))
#> DataFrame with 6 rows and 10 columns
#>                                                                original_name    database   size_db size_universe                           contributions
#>                                                                  <character> <character> <integer>     <integer>                                  <list>
#> kegg.hsa04080 Neuroactive ligand-receptor interaction hsa04080 Neuroactive..        kegg       370           195  0.00508377,-0.04229206,-0.02369940,...
#> kegg.hsa04142 Lysosome biogenesis                     hsa04142 Lysosome bi..        kegg       250           237    -0.0183027,-0.0672571,-0.0180544,...
#> kegg.hsa04151 PI3K-Akt signaling pathway              hsa04151 PI3K-Akt si..        kegg       362           300  0.00888585,-0.00314775,-0.00335746,...
#> kegg.hsa04310 Wnt signaling pathway                   hsa04310 Wnt signali..        kegg       178           152  0.04872311,-0.00818062,-0.03842735,...
#> kegg.hsa04517 IgSF CAM signaling                      hsa04517 IgSF CAM si..        kegg       299           270 -0.11821237,-0.00722345, 0.05590336,...
#> kegg.hsa04519 Cadherin signaling                      hsa04519 Cadherin si..        kegg       404           291    -0.1348806,-0.0715228,-0.0366950,...
#>                                                           baseMean log2FoldChange        pval        padj        sd
#>                                                          <numeric>      <numeric>   <numeric>   <numeric> <numeric>
#> kegg.hsa04080 Neuroactive ligand-receptor interaction -6.97854e-17       -9.07018 8.36299e-25 1.28026e-22   4.10314
#> kegg.hsa04142 Lysosome biogenesis                     -2.09753e-16        8.50997 2.14148e-16 2.00996e-15   4.03683
#> kegg.hsa04151 PI3K-Akt signaling pathway              -5.59870e-16       -9.86022 3.46713e-20 9.35703e-19   4.54484
#> kegg.hsa04310 Wnt signaling pathway                    1.47699e-17       -8.13649 4.50879e-21 1.68208e-19   3.73294
#> kegg.hsa04517 IgSF CAM signaling                       2.12528e-16       -9.26569 7.09007e-19 1.23437e-17   4.30560
#> kegg.hsa04519 Cadherin signaling                       6.32827e-16       -8.63763 2.95062e-14 1.80161e-13   4.20368
```

We can visualize the result as a heatmap with the contribution of each
gene to the score.

``` r

GSDS.Heatmap(resGSDS.sel, colData = sampleAnnot["culture_media"])
```

<img src="ReadMe_files/figure-gfm/heatmapGSDS-1.png" alt="" width="100%" style="display: block; margin: auto;" />


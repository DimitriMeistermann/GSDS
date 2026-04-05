## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(7,7),
  fig.align = "center",
  cache = TRUE,
  warning=FALSE
)

## ----FCS----------------------------------------------------------------------
fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
resEnrich<-enrich.fcs(fcsScore,DBsets = geneSetDB,returnGenes = TRUE) #return genes of the gene set in the dataframe
head(resEnrich[order(resEnrich$pval,decreasing = FALSE),])

# The result can be visualized as a volcano plot

volcanoPlot(resEnrich,effectSizeCol = "NES", adjPvalCol = "padj",labelCol = resEnrich$term)


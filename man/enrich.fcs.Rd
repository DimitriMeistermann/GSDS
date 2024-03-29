% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeEnrichment.R
\name{enrich.fcs}
\alias{enrich.fcs}
\title{Functional class scoring enrichment (fgsea algorithm)}
\usage{
enrich.fcs(
  x,
  idGeneDF = NULL,
  database = c("kegg", "reactom", "goBP", "goCC", "goMF"),
  maxSize = 500,
  minSize = 2,
  customAnnot = NULL,
  returnGenes = FALSE,
  keggDisease = FALSE,
  species = "Human",
  db_terms = NULL,
  speciesData = NULL,
  ...
)
}
\arguments{
\item{x}{vector or dataframe/matrix of one column. Values are used to classify the genes, example, it can be Log2(Fold-Change0). Genes are contained in the names/rownames of the vector/dataframe. Example of valid x: x<-rnorm(n = 4); names(x)<-c("GATA2","SOX17","KLF4","POU5F1")}

\item{idGeneDF}{Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL \code{species} and \code{speciesData} arguments wont be used.}

\item{database}{Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom}

\item{maxSize}{maximum number of gene in each term}

\item{minSize}{Minimum number of gene in each term}

\item{customAnnot}{custom annotation database, as a list of terms, each element contain a vector of gene symbols.}

\item{returnGenes}{return genes that were the most important for the enrichment of term}

\item{keggDisease}{Logical. Retain kegg disease term ?}

\item{species}{Character. Shortname of the species as described in \code{data("bods")}.}

\item{db_terms}{A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
If this argument is not NULL, no additional database are downloaded.}

\item{speciesData}{object returned by \code{getSpeciesData}. If not NULL \code{species} argument wont be used.}

\item{...}{Additional parameters that are passed to fgsea}
}
\value{
A data frame with the following columns:
\itemize{
\item pathway: name of the pathway/term
\item pval: an enrichment p-value
\item padj: a BH-adjusted p-value
\item log2err: the expected error for the standard deviation of the P-value logarithm
\item ES: enrichment score, same as in Broad GSEA implementation
\item NES: enrichment score normalized to mean enrichment of random samples of the same size
\item size: number of gene in the term after removing genes not present
\item genes (if \code{returnGenes}). Vector of genes of the term.
}
}
\description{
Functional class scoring enrichment (fgsea algorithm)
}
\examples{
data("DEgenesPrime_Naive",package = "oob")
fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),
    DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
resEnrich<-enrich.fcs(fcsScore,database = "kegg",species = "Human")
head(resEnrich)
}

activeScorePCAlist<-function (exprMatrix, geneList, scale = FALSE, center = TRUE)
{
    res <- lapply(geneList, function(genesOfEl) activScorePCA(exprMatrix,
                                                              genesOfEl, returnContribution = TRUE, scale = scale,
                                                              center = center))
    list(activScoreMat = sapply(res, function(x) x$activScore),
         contributionList = lapply(res, function(x) x$contribution))
}



#' Compute the activation score of gene sets from an expression matrix.
#'
#' @description Perform a PCA for each gene set, from the matrix of [ genes from gene set × all samples ]. Return the first PCs as activations scores of the gene sets.
#'
#' @param expressionMatrix An expression matrix (normalized log2(x+1) counts). Genes as rows and sample as columns. If `db_terms` is not given, must be named by gene symbols.
#' @param corrIdGenes Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param scaleScores Logical. Divide expression of gene by its standard deviation before doing the PCA.
#' @param centerScores Logical. Subtract mean to gene expression before doing the PCA.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom.
#' @param maxSize Maximum number of gene in each term.
#' @param minSize Minimum number of gene in each term.
#' @param customAnnot Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param db_terms A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#'
#' @return A list where each element is a database of gene set given as input.
#' For each database, contain a list of activation score, with gene sets as rows and samples as columns ;
#' and the list of contribution (or weight) to activation score of each gene per gene set.
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' keggDB<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,db_terms = keggDB)
#' #same as
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,database = "kegg")
computeActivationScore<-function(expressionMatrix,corrIdGenes=NULL,scaleScores=FALSE,centerScores=TRUE,
                                 database=c("kegg","reactom","goBP","goCC","goMF"),
                                 maxSize=500,minSize=2,customAnnot=NULL,
                                 keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL){

    if(!class(expressionMatrix)[1]%in%c("data.frame","matrix")) stop("expressionMatrix should be a matrix or a dataframe")
    if(class(rownames(expressionMatrix))!="character") stop("rows of expression matrix should be named with genes symbol")
    if(is.null(db_terms)){
        db_terms<-getDBterms(geneSym=rownames(expressionMatrix), corrIdGenes=corrIdGenes,database=database,
                             customAnnot=customAnnot,keggDisease=keggDisease,species=species,returnGenesSymbol = TRUE)
    }

    if(length(db_terms)==0) stop("Error, no term in any database was found")

    lapply(db_terms,function(database){
        database<-lapply(database,function(genesOfTerm) intersect(genesOfTerm,rownames(expressionMatrix)))
        nGenePerTerm<-sapply(database,length)
        database<-database[nGenePerTerm>minSize & nGenePerTerm<maxSize]

        return(activeScorePCAlist(exprMatrix, database, scale = scaleScores))
    })
}



#' Compute a "interest score" for a set of genes. Useful for GSEA.
#'
#' @param genes A vector of gene names
#' @param pvalues A vector of numeric corresponding to p-values.
#' @param logFoldChanges A vector of numeric corresponding to Log2(Fold-change)
#' @param logPval Logical. Compute the p-value score as `-log10(pval)` instead of `1-pval`.
#'
#' @return A vector of numeric corresponding to interest scores, named by genes.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
fcsScoreDEgenes<-function(genes,pvalues,logFoldChanges,logPval=FALSE){
    if( sum(length(genes) == c(length(pvalues),length(logFoldChanges))) <2) stop("genes, pvalues and logPval should have the same length")
    if(logPval){
        pvalScore<- -log10(pvalues)
    }else{
        pvalScore <- 1-pvalues
    }
    pvalScore[logFoldChanges<0]<- -pvalScore[logFoldChanges<0]
    names(pvalScore)<-genes
    pvalScore
}

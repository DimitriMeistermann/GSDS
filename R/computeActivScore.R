

#' Compute the activation score of a gene set from 1st component of its PCA
#'
#' @param data A matrix of numeric with rows as features (in the RNA-Seq
#'   context, log counts). Can also be a SingleCellExperiment object.
#' @param genes A character vector. The gene set where the activation score has
#'   to be computed. Must be a subset of `data` row names.
#' @param transpose Logical. If TRUE, `data` is transposed with `t()` before
#'   computing the PCA.
#' @param scale  Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#' @param returnContribution Logical. Return list with activation score and
#'   contribution of genes to the activation score.
#' @param sce_assay Integer
#'   or character, if `data` is a `SingleCellExperiment` object, the assay name
#'   to use.
#' @return A vector of numeric corresponding to activation scores, named by genes. If `returnContribution` return a list with activation scores and contributions of genes.
#' @export
#'
#' @examples
#' data("bulkLogCounts", package = "oob")
#' keggData<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' geneSet<-keggData$kegg$`hsa00190 Oxidative phosphorylation`
#' geneSet<-intersect(geneSet,rownames(bulkLogCounts))
#' activScorePC1(bulkLogCounts,genes = geneSet)
#' activScorePC1(bulkLogCounts,genes = geneSet,returnContribution = TRUE)
activScorePC1 <- function (data,
													 genes,
													 transpose = TRUE,
													 scale = FALSE,
													 center = TRUE,
													 returnContribution = FALSE,
													 sce_assay = "logcounts")
{
	if (inherits(data, "SingleCellExperiment")) {
		sce_obj <- data
		data <- assay(sce_obj, sce_assay)
	}
	if (transpose)
		data <- t(data)
	if (sum(!genes %in% colnames(data)) > 0)
		stop("genes should be a subset of data row names")
	pca <- fastPCA(
		data[, genes],
		center = center,
		scale = scale,
		nPC = 1,
		transpose = FALSE
	)
	activScore <- pca$x[, 1]
	contribution <- pca$rotation[, 1]
	if (cor(rowMeans(data[, genes]), activScore) < 0) {
		activScore <- -activScore
		contribution <- -contribution
	}
	if (returnContribution) {
		list(activScore = activScore, contribution = contribution)
	}
	else {
		activScore
	}
}

#' Compute activation score for a list of gene sets.
#'
#' @param data A matrix of expression values, genes as rows and samples as columns. Can also be a SingleCellExperiment object.
#' @param geneList A list of gene sets, each element is a vector of genes that exist in `data` row names.
#' @param transpose Logical. If TRUE, `data` is transposed with `t()` before computing the PCA.
#' @param scale Logical. Divide expression of gene by its standard deviation before doing the PCA.
#' @param center Logical. Subtract mean to gene expression before doing the PCA.
#' @param transpose Logical. If TRUE, `data` is transposed with `t()` before computing the PCA.
#' @param sce_assay Integer
#'   or character, if `data` is a `SingleCellExperiment` object, the assay name
#'   to use.
#' @return A list with two elements:
#' - activScoreMat: A matrix of activation score, with gene sets as rows and samples as columns.
#' - contributionList: A list of contribution (or weight) to activation score of each gene per gene set.
#' @export
#'
#' @examples
#' data("bulkLogCounts", package = "oob")
#' keggDB<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' geneSetActivScore<-activScorePC1list(bulkLogCounts,geneList = keggDB[[1]])
activScorePC1list <-
    function (data,
              geneList,
              transpose = TRUE,
              scale = FALSE,
              center = TRUE,
    					sce_assay = "logcounts")
    {
	    	if (inherits(data, "SingleCellExperiment")) {
	    		sce_obj <- data
	    		data <- assay(sce_obj, sce_assay)
	    	}
        if (!transpose) data <- t(data)
        sd_pos <- apply(data,2,sd) > 0
				if(sum(!sd_pos) > 0) warning("Some features have sd = 0, they will be deleted from the analysis")
        data <- data[,sd_pos]
        res <-
            lapply(geneList, function(genesOfEl)
                activScorePC1(
                    data,
                    genesOfEl,
                    returnContribution = TRUE,
                    scale = scale,
                    center = center
                ))

        for (genesOfEl in geneList)
            activScorePC1(
                data,
                genesOfEl,
                returnContribution = TRUE,
                scale = scale,
                center = center
            )
        list(
            activScoreMat = sapply(res, function(x)
                x$activScore),
            contributionList = lapply(res, function(x)
                x$contribution)
        )
    }



#' Compute the activation score of gene sets from an expression matrix.
#'
#' @description Perform a PCA for each gene set, from the matrix of *genes from
#'   gene set × all samples*. Return the first PCs as activation scores of the
#'   gene sets.
#'
#' @param data An expression matrix (normalized log2(x+1) counts). Genes as rows
#'   and sample as columns by default. If `db_terms` is not given, must be named by gene
#'   symbols. Can also be a SingleCellExperiment object.
#' @param transpose Logical. If `TRUE`, `data` is transposed with `t()` before
#'   computing the PCA.
#' @param idGeneDF Dataframe of gene ID correspondence where each column is a
#'   gene ID type. If not NULL `species` and `speciesData` arguments wont be
#'   used.
#' @param scaleScores Logical. Divide expression of gene by its standard
#'   deviation before doing the PCA.
#' @param centerScores Logical. Subtract mean to gene expression before doing
#'   the PCA.
#' @param database Which annotation database ? valid: database: kegg reactom
#'   goBP goCC goMF custom.
#' @param maxSize Maximum number of gene in each term.
#' @param minSize Minimum number of gene in each term.
#' @param customAnnot Custom annotation database, as a list of terms, each
#'   element contain a vector of gene symbols.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in
#'   `data("bods")`.
#' @param db_terms A list or NULL. A named list were each element is a database.
#'   Inside each database, a list terms, named by the term and containing gene
#'   vectors as gene symbols.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species`
#'   argument wont be used.
#' @return A list where each element is a database of gene set given as input.
#'   For each database, contain a list of activation score, with gene sets as
#'   rows and samples as columns ; and the list of contribution (or weight) to
#'   activation score of each gene per gene set.
#' @param sce_assay Integer or character, if `data` is a `SingleCellExperiment`
#'   object, the assay name to use.
#' @export
#'
#' @examples
#' data("bulkLogCounts", package = "oob")
#' keggDB<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,db_terms = keggDB)
#' #same as
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,database = "kegg")
computeActivationScore <-
    function(data,
             transpose = TRUE,
             idGeneDF = NULL,
             scaleScores = FALSE,
             centerScores = TRUE,
             database = c("kegg", "reactom", "goBP", "goCC", "goMF"),
             maxSize = 500,
             minSize = 2,
             customAnnot = NULL,
             keggDisease = FALSE,
             species = "Human",
             db_terms = NULL,
             speciesData = NULL,
    				 sce_assay = "logcounts") {
	    	if (inherits(data, "SingleCellExperiment")) {
	    		sce_obj <- data
	    		data <- assay(sce_obj, sce_assay)
	    	}
        if ( !is.data.frame(data) & !is.matrix(data))
            stop("data should be a matrix or a dataframe")
        if(!transpose)
            data<-t(data)
        if (! is.character(rownames(data)))
            stop("rows of expression matrix should be named with genes symbol")

        valid <- apply(data,1,var) > 0
				if(sum(!valid)>0) warning(length(valid)," feature(s) with 0 variance, they will be removed from the input")

        if (is.null(db_terms)) {
            db_terms <-
                getDBterms(
                    geneSym = rownames(data),
                    idGeneDF = idGeneDF,
                    database = database,
                    customAnnot = customAnnot,
                    keggDisease = keggDisease,
                    species = species,
                    returnGenesSymbol = TRUE
                )
        }
        if (length(db_terms) == 0)
            stop("Error, no term in any database was found")

        lapply(db_terms, function(database) {
            database <-
                lapply(database, function(genesOfTerm)
                    intersect(genesOfTerm, rownames(data)))
            nGenePerTerm <- sapply(database, length)
            database <-
                database[nGenePerTerm > minSize & nGenePerTerm < maxSize]

            return(activScorePC1list(data, database, scale = scaleScores))
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
#' data("DEgenesPrime_Naive", package = "oob")
#' fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),
#'     DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
fcsScoreDEgenes <-
    function(genes, pvalues, logFoldChanges, logPval = FALSE) {
        if (sum(length(genes) == c(length(pvalues), length(logFoldChanges))) <
            2)
            stop("genes, pvalues and logPval should have the same length")
        if (logPval) {
            pvalScore <- -log10(pvalues)
        } else{
            pvalScore <- 1 - pvalues
        }
        pvalScore[logFoldChanges < 0] <- -pvalScore[logFoldChanges < 0]
        names(pvalScore) <- genes
        pvalScore
    }


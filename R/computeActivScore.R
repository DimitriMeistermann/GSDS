

#' Compute the activation score of a gene set from 1st component of its PCA
#'
#' @param data A matrix of numeric with rows as features (in the RNA-Seq
#'   context, log counts). Can also be a SummarizedExperiment object.
#' @param genes A character vector. The gene set where the activation score has
#'   to be computed. Must be a subset of `data` row names.
#' @param transpose Logical. If TRUE, `data` is transposed with `t()` before
#'   computing the PCA.
#' @param scale  Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#' @param returnContribution Logical. Return list with activation score and
#'   contribution of genes to the activation score.
#' @param se_data_assay Integer
#'   or character, if `data` is a `SummarizedExperiment` object, the assay name
#'   to use.
#' @return A vector of numeric corresponding to activation scores, named by genes. If `returnContribution` return a list with activation scores and contributions of genes.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("keggHuman")
#' geneSet<-keggHuman$kegg$`hsa00190 Oxidative phosphorylation`
#' geneSet<-intersect(geneSet,rownames(bulkLogCounts))
#' activScorePC1(bulkLogCounts,genes = geneSet)
#' activScorePC1(bulkLogCounts,genes = geneSet,returnContribution = TRUE)
activScorePC1 <- function (data,
													 genes,
													 transpose = TRUE,
													 scale = FALSE,
													 center = TRUE,
													 returnContribution = FALSE,
													 se_data_assay = "logcounts")
{
	if (inherits(data, "SummarizedExperiment")) {
		sce_obj <- data
		data <- assay(sce_obj, se_data_assay)
	}
	if (transpose)
		data <- Matrix::t(data)
	if (sum(!genes %in% colnames(data)) > 0)
		stop("genes should be a subset of data row names")
	pca <- fastPCA(
		data[, genes, drop = FALSE],
		center = center,
		scale = scale,
		nPC = 1,
		transpose = FALSE
	)
	activScore <- pca$x[, 1]
	contribution <- pca$rotation[, 1]
	correlation_val <- cor(Matrix::rowMeans(data[, genes, drop = FALSE]), activScore)
	if (!is.na(correlation_val) && correlation_val < 0) {
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
#' @inheritParams computeActivationScore
#' @param data A matrix of expression values, genes as rows and samples as columns. Can also be a SummarizedExperiment object.
#' @param geneList A list of gene sets, each element is a vector of genes that exist in `data` row names.
#' @param transpose Logical. If TRUE, `data` is transposed with `t()` before computing the PCA.
#' @param scale Logical. Divide expression of gene by its standard deviation before doing the PCA.
#' @param center Logical. Subtract mean to gene expression before doing the PCA.
#' @param transpose Logical. If TRUE, `data` is transposed with `t()` before computing the PCA.
#' @return A list with two elements:
#' - activScoreMat: A matrix of activation score, with gene sets as rows and samples as columns.
#' - contributionList: A list of contribution (or weight) to activation score of each gene per gene set.
#'
#' @examples
#' data("bulkLogCounts")
#' data("keggHuman")
#' geneList <-  keggHuman$kegg[seq.int(5)]
#' geneList <- lapply(geneList, function(x) intersect(x,rownames(bulkLogCounts)))
#' geneSetActivScore<-activScorePC1list(bulkLogCounts,geneList = geneList)
activScorePC1list <- function(data,
															geneList,
															transpose = TRUE,
															scale = FALSE,
															center = TRUE,
															BPPARAM = BiocParallel::SerialParam()) {

	if (transpose) data <- Matrix::t(data)

	res <- BiocParallel::bplapply(geneList, function(genesOfEl) {
		# Ensure only relevant genes are passed to the inner function to save memory
		suppressWarnings(activScorePC1(
			data = data[,genesOfEl],
			genes = genesOfEl,
			returnContribution = TRUE,
			scale = scale,
			center = center,
			transpose = FALSE
		))
	}, BPPARAM = BPPARAM)

	return(list(
		activScoreMat = Matrix::t(vapply(res, function(x) x$activScore, numeric(nrow(data)))),
		contributionList = lapply(res, function(x) x$contribution)
	))
}


#' Compute the activation score of gene sets from an expression matrix.
#'
#' @description Perform a PCA for each gene set, from the matrix of *genes from
#'   gene set × all samples*. Return the first PCs as activation scores of the
#'   gene sets.
#'
#' @param data An expression matrix (normalized log2(x+1) counts). Genes as rows
#'   and sample as columns by default. If `db_terms` is not given, must be named
#'   by gene symbols. Can also be a SummarizedExperiment object.
#' @param DBsets A list or NULL. A named list were each element is a database.
#'   Inside each database, a list terms, named by the term and containing gene
#'   vectors as gene symbols.
#' @param transpose Logical. If `TRUE`, `data` is transposed with `t()` before
#'   computing the PCA.
#' @param scaleScores Logical. Divide expression of gene by its standard
#'   deviation before doing the PCA.
#' @param maxSize Maximum number of gene in each term.
#' @param minSize Minimum number of gene in each term.
#' @param se_data_assay Integer
#'   or character, if `data` inherit a `SummarizedExperiment` object, the assay name where the input matrix is.
#' @param BPPARAM Parallelization parameter as used in BiocParallel::bpparam().
#' @param returnSummarizedExperiment do not keep the results listed by database but merge them as a single activation SummarizedExperiment object
#'
#' @return if `returnSummarizedExperiment = TRUE`
#' A SummarizedExperiment object, with the activation score as the assay, and gene set metadata stored in the `rowData`.
#' The gene set metadata contains size_db (number of gene in the term originally present in the given database),
#' size_universe (number of gene in the term after removing genes not present).
#' if `returnSummarizedExperiment = FALSE`
#' A list by database containing each a list with the matrix of activation score, with gene sets as
#'   rowa and samples as columns ; and the list of contribution (or weight) to
#'   activation score of each gene per gene set, these two elements are superseded by the list of database given as input.
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("keggHuman")
#' geneSetActivScore<-computeActivationScore(bulkLogCounts, keggHuman)
computeActivationScore <- function(data,
																	 DBsets,
																	 transpose = TRUE,
																	 scaleScores = FALSE,
																	 returnSummarizedExperiment = TRUE,
																	 maxSize = 500,
																	 minSize = 2,
																	 se_data_assay = "logcounts",
																	 BPPARAM = BiocParallel::SerialParam()) {

	if(!checkDB(DBsets)) stop("DBsets is empty or not valid. Pease check the databases contains a list of characters")
	if (inherits(data, "SummarizedExperiment")) {
		data <- SummarizedExperiment::assay(data, se_data_assay)
	}

	if (transpose) data <- Matrix::t(data) # Ensure rows = Genes, cols = Cells
	if(is.data.frame(data)) data <- as.matrix(data)

	rvars <- if(inherits(data, "dgCMatrix")) {
		sparseMatrixStats::colVars(data)
	} else {
		matrixStats::colVars(data)
	}

	valid <- rvars > 0
	if (any(!valid)) {
		warning(sum(!valid), " feature(s) with 0 variance removed.")
		data <- data[, valid, drop = FALSE]
	}
	valid_genes <- colnames(data)
	res <- lapply(DBsets, function(database) {
		size_db <- lengths(database)
		database <- lapply(database, function(genes) genes[genes %in% valid_genes])
		size_universe <- lengths(database)
		selDB <- size_universe >= minSize & size_universe <= maxSize
		activScore <- activScorePC1list(data, database[selDB], transpose = FALSE, scale = scaleScores,BPPARAM = BPPARAM)
		attr(activScore,"size_db") <- size_db[selDB]
		attr(activScore,"size_universe") <- size_universe[selDB]
		return(activScore)
	})
	res <- res[!sapply(res, is.null)] # Remove empty results
	if (!returnSummarizedExperiment) return(res)
	db_lengths <- sapply(res, function(x) nrow(x$activScoreMat))
	db_names <- rep(names(res), times = db_lengths)
	original_names <- unlist(lapply(res, function(x) rownames(x$activScoreMat)), use.names = FALSE)
	renamed_terms <- paste0(db_names, ".", original_names)

	# Combine Matrices (Cells x Pathways)
	combined_scores  <- do.call(rbind, lapply(res, `[[`, "activScoreMat"))
	# Combine Contributions (List of vectors)
	combined_contrib <- do.call(c, lapply(res, `[[`, "contributionList"))

	se <- SummarizedExperiment::SummarizedExperiment(
		assays  = list(activScoreMat = combined_scores),
		rowData = S4Vectors::DataFrame(
			original_name = original_names,
			database      = db_names,
			size_db       = lapply(res,attr,"size_db") |> unlist(),
			size_universe = lapply(res,attr,"size_universe") |> unlist(),
			contributions = I(combined_contrib)
		)
	)
	colnames(se) <- rownames(data)
	rownames(se) <- renamed_terms

	return(se)
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


fastPCA<-function(data,
				 transpose = TRUE,
				 scale = FALSE,
				 center = TRUE,
				 nPC = min(ncol(data) - 1, nrow(data) - 1, 30),
				 weight.by.var = TRUE,
				 ...) {


	if(!is.matrix(data)) data <- as.matrix(data)
	if (transpose)
		data <- t(data)
	means <- 0
	sdeviations <- 1
	if (center | scale)
		data <- scale(data, scale = scale, center = center)
	if (center)
		means <- attr(data, "scaled:center")
	if (scale)
		sdeviations <- attr(data, "scaled:scale")

	resacp <- list()
	resacp$n.obs <- dim(data)[1]

	resacp$scale <- scale
	resacp$center <- center
	resacp$transform <-
		list(sdeviations = sdeviations, means = means)

	irlbaResults <- irlba::irlba(A = data, nv = nPC, ...)
	rotation <- irlbaResults$v
	resacp$sdev <- irlbaResults$d / sqrt(max(1, nrow(data) - 1))
	if (weight.by.var) {
		if (nPC > 1) {
			reducedSpace <- irlbaResults$u %*% diag(irlbaResults$d)
		} else {
			reducedSpace <- irlbaResults$u %*% irlbaResults$d
		}
	} else {
		reducedSpace <- irlbaResults$u
	}
	rownames(rotation) <- colnames(data)
	colnames(rotation) <- paste0("PC", seq_len(nPC))
	rownames(reducedSpace) <- rownames(data)
	colnames(reducedSpace) <- colnames(rotation)
	resacp$x <- reducedSpace
	resacp$rotation <- rotation
	resacp$propExplVar <- resacp$sdev ^ 2 / sum(resacp$sdev ^ 2)
	resacp$isFastPCA <- TRUE
	return(resacp)
}

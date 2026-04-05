#' Functional class scoring enrichment (fgsea algorithm)
#'
#' @inheritParams computeActivationScore
#' @param x vector or dataframe/matrix of one column. Values are used to classify the genes, example, it can be Log2(Fold-Change0). Genes are contained in the names/rownames of the vector/dataframe. Example of valid x: x<-rnorm(n = 4); names(x)<-c("GATA2","SOX17","KLF4","POU5F1")
#' @param returnGenes  return genes that were the most important for the enrichment of term
#' @param ... Additional parameters that are passed to fgsea
#'
#' @return
#' A data frame with the following columns:
#' - term: name of the term/gene set
#' - database: origin of the gene set
#' - size: number of gene in the term after removing genes not present.
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#' - log2err: the expected error for the standard deviation of the P-value logarithm
#' - ES: enrichment score, same as in Broad GSEA implementation
#' - NES: enrichment score normalized to mean enrichment of random samples of the same size
#' - genes (if `returnGenes`). Vector of genes of the term.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' data("keggHuman")
#' fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),
#'     DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
#' resEnrich<-enrich.fcs(fcsScore,DBsets =  keggHuman)
#' head(resEnrich)
enrich.fcs<-function(x, DBsets,
                     maxSize=500,minSize=2,returnGenes=FALSE,
										 BPPARAM = BiocParallel::SerialParam(),...){
    if(is.data.frame(x) | is.matrix(x)){
        if(is.null(rownames(x))) stop("If input format is a matrix/dataframe, row should be named with gene symbols")
        tempx<-x
        x<-tempx[,1]
        names(x)<-rownames(tempx)
    }

    if(is.null(names(x))) stop("Values should be named with gene symbols")
    if(!is.numeric(x)) stop("Values must be numeric")
    if(!checkDB(DBsets)) stop("DBsets is empty or not valid. Pease check the databases contains a list of characters")
    res<-list()
    for(db in names(DBsets)){
    	resGSEA <- suppressWarnings(fgsea::fgseaMultilevel (
    		DBsets[[db]],
    		x , BPPARAM = BPPARAM,
    		minSize = minSize,
    		maxSize = maxSize,
    		eps = 0,
    		...
    	))
        resGSEA$leadingEdge<-NULL
        colnames(resGSEA)[1]<-"term"
        res[[db]] <-
            data.frame(resGSEA[1],
                       "database" = db,
                       "size" = sapply(DBsets[[db]][resGSEA$term], length),
                       resGSEA[seq(2, 6)])
        if(returnGenes) res[[db]]$genes <- DBsets[[db]][res[[db]]$term]
    }
    res<-do.call("rbind", res)
    res$padj<-p.adjust(res$pval,method = "BH")
    return(data.frame(res))
}

#to do, make it compatible with multiple vectors for venn diagram like enrichment

#' Over Representation Analysis from Hypergeometric test
#'
#' @inheritParams computeActivationScore
#' @inheritParams enrich.fcs
#' @param x vector or dataframe/matrix of one column.
#' Values are booleans and say if gene is from the list of interest or not.
#' Genes are contained in the names/rownames of the vector/dataframe.
#' Example of valid x: x<-c(TRUE,TRUE,FALSE,FALSE); names(x)<-c("GATA2","SOX17","KLF4","POU5F1").
#' In this case, GATA2, SOX17, KLF4, POU5F1 are the universe of gene and GATA2 and SOX17 are the genes of interest
#'
#' @return
#' A dataframe with the following columns:
#' - term: name of the term/gene set
#' - database: origin of the gene set
#' - size_db: number of gene in the term originally present in the given database.
#' - size_universe: number of gene in the term after removing genes not present.
#' - obsOverlap: number of gene in common between set of interest and term
#' - expOverlap: expected number of gene in common between set of interest and term
#' - OEdeviation: \eqn{\frac{obs - exp}{\sqrt{Universe}}}, normalized deviation from expected value
#' - log2OR: log2 of the Odds Ratio
#' - pval: p-value based on a Hypergeometric test
#' - padj: a BH-adjusted p-value
#' - genes (if `returnGenes`). Vector of genes of the term.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' data("keggHuman")
#' vectorIsDE<-DEgenesPrime_Naive$isDE!="NONE";names(vectorIsDE)<-rownames(DEgenesPrime_Naive)
#' resEnrich<-enrich.ora(vectorIsDE, DBsets = keggHuman)
#' head(resEnrich)
enrich.ora <-
	function(x, DBsets,
					 minSize = 2,
					 maxSize = 500,
					 returnGenes = FALSE) {
		if (is.data.frame(x) | is.matrix(x)) {
			tempx <- x
			x <- tempx[, 1]
			names(x) <- rownames(tempx)
		}

		if (!is.logical(x))
			stop("Values must be logical (TRUE or FALSE)")
		if (!is.character(names(x)))
			stop("Values must be named with gene IDs")

		if(!checkDB(DBsets)) stop("DBsets is empty or not valid. Pease check the databases contains a list of characters")

		nInterest <- sum(x)
		nuniverse <- length(x)

		results <- list()
		for (db in names(DBsets)) {
			terms <- DBsets[[db]]

			nGeneBytermTotal <- sapply(terms, length)
			terms <- sapply(terms, function(term) {
				return(intersect(term, names(x)))
			})
			nGeneByterm<-sapply(terms, length)

			goodterm <- nGeneByterm >= minSize & nGeneByterm <= maxSize
			terms <- terms[goodterm]
			nGeneByterm <- nGeneByterm[goodterm]
			nGeneBytermTotal <- nGeneBytermTotal[goodterm]

			nGeneOfInterestByterm <- sapply(terms, function(term) {
				return(sum(x[term]))
			})


			results[[db]] <- data.frame(row.names = names(terms))
			results[[db]]$term <- names(terms)
			results[[db]]$database <- db
			results[[db]]$size_db <- nGeneBytermTotal

			parameterList4Enrich <-
				vector(mode = "list", length = length(terms))
			for (i in seq_along(terms)) {
				parameterList4Enrich[[i]] <-
					list(
						intersectionSize = nGeneOfInterestByterm[i],
						setSizes = c(nInterest, nGeneByterm[i]),
						universeSize = nuniverse
					)
			}

			resEnrich <-
				lapply(parameterList4Enrich, function(params)
					do.call("getORAstats", params))
			results[[db]]$size_universe <- nGeneByterm

			results[[db]]$obsOverlap <- nGeneOfInterestByterm
			results[[db]]$expectOverlap <-
				sapply(resEnrich, function(x)
					x$expected)
			results[[db]]$OEdeviation <-
				sapply(resEnrich, function(x)
					x$OEdeviation)
			results[[db]]$log2OR <- sapply(resEnrich, function(x)
				x$log2OR)
			results[[db]]$pval <- sapply(resEnrich, function(x)
				x$pval)
			results[[db]]$padj <- 0

			if (returnGenes) {
				results[[db]]$genes <- DBsets[[db]]
			}
		}
		results <- do.call("rbind", results)
		results$padj <- p.adjust(results$pval, method = "BH")
		return(results)
	}

#' Gene Set Differential Scoring (GSDS)
#' @inheritParams computeActivationScore
#' @param geneSetActivScore A `SummarizedExperiment` as returned by `computeActivationScore`
#' @param colData An annotation dataframe. Each column is a feature, each row a sample. Same number of samples than in `data`.
#' @param contrast A vector of 3 character.
#' 1. Name of the experimental variable that have to be used for differential activation. Must be a column name of `colData`.
#' 2. Condition considered as the reference.
#' 3. Condition considered as the target group.
#' @param formula Formula for the linear modelling, if NULL formula is group. Example : `~ Condition + Batch` if condition and
#'   Batch are colnames of colData
#' @param return_seobject return a `SummarizedExperiment` object, with rowData annotated with GSDS results instead of the result dataframe.
#'   to use.
#' @return
#' A dataframe name by gene set IDs, with the following columns:
#' - original_name: name of the term/gene set as it is in the provided database
#' - database: origin of the gene set
#' - size_db: number of gene in the gene set.
#' - size_universe:number of gene in the gene set that are presesent in the data (e.g. used for computing the activation score)
#' - baseMean: mean of activation score in the gene set
#' - sd: standard deviation  in the gene set
#' - log2FoldChange: Log(Log Fold Change) of activation score between the two tested groups.
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("keggHuman")
#'
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,DBsets = keggHuman)
#' resGSDS<-GSDS(geneSetActivScore = geneSetActivScore,colData = sampleAnnot,
#'     contrast = c("culture_media","T2iLGO","KSR+FGF2"))
#' head(rowData(resGSDS))
#' # correct for the line effect
#' resGSDS<-GSDS(geneSetActivScore = geneSetActivScore,colData = sampleAnnot,
#'     contrast = c("culture_media","T2iLGO","KSR+FGF2"), formula = ~ culture_media + line)
#' head(rowData(resGSDS))
#'
#' #just a data frame
#' head(GSDS(geneSetActivScore = geneSetActivScore,colData = sampleAnnot,
#'     contrast = c("culture_media","T2iLGO","KSR+FGF2"),DBsets =  keggHuman,
#'     return_seobject = FALSE))
#' # or
#' resGSDS<-GSDS(data = bulkLogCounts,colData = sampleAnnot,
#'     contrast = c("culture_media","T2iLGO","KSR+FGF2"),DBsets = keggHuman)
#'
GSDS <-
    function(geneSetActivScore=NULL, DBsets=NULL,
             data = NULL,
    				 colData,
             contrast ,
    				 formula = NULL,
             maxSize = 500,
             minSize = 2,
             se_data_assay = "logcounts",
    				 return_seobject = TRUE,
    				 BPPARAM = BiocParallel::SerialParam()) {

        if (is.null(geneSetActivScore)){
        	if(is.null(data)) stop("At least data or geneSetActivScore must be given")
        	if(is.null(DBsets)) stop("If geneSetActivScore is not given, DBsets must be set")
        	geneSetActivScore <-
        		computeActivationScore(data = data, DBsets = DBsets)
        } else{
        	invisible(checkGSDSobj(geneSetActivScore,"geneSetActivScore"))
        }

        if (!setequal(rownames(colData),
                      colnames(geneSetActivScore))) {
            stop("colData and data must have the same samples names, in same order")
        }
    	activScoreMat <- assay(geneSetActivScore,"activScoreMat")
    	if(return_seobject){
    		GSDSdt <- data.frame(multiLinearModel(activScoreMat, colData, contrast = contrast, formula=formula),
    												 sd = apply(activScoreMat, 1, sd))
    		rowData(geneSetActivScore) <- cbind(rowData(geneSetActivScore), GSDSdt)
    		return(geneSetActivScore)
    	}else{
    		metaData <- rowData(geneSetActivScore)
    		metaData$contributions <- NULL
    		return(data.frame(metaData, multiLinearModel(activScoreMat, colData, contrast = contrast),
    										sd = apply(activScoreMat, 1, sd)))
    	}
    }

multiLinearModel <- function(exprData,
														 colData,
														 contrast,
														 formula = NULL) {

	target_col <- contrast[1]
	valid_samples <- colData[[target_col]] %in% contrast[2:3]

	sub_expr <- exprData[, valid_samples, drop = FALSE]
	sub_meta <- colData[valid_samples, , drop = FALSE]


	sub_meta[[target_col]] <- factor(sub_meta[[target_col]], levels = contrast[3:2])

	if (is.null(formula)) {
		# Default to simple group comparison
		test_formula <- as.formula(paste("y ~", target_col))
	} else {
		form_char <- as.character(formula)
		rhs <- form_char[length(form_char)]
		test_formula <- as.formula(paste("y ~", rhs))
	}

	resList <- lapply(seq_len(nrow(sub_expr)), function(i) {

		regTab <- sub_meta
		regTab$y <- sub_expr[i, ]

		# Defensive Check: Try to fit the model, catch errors if cases are 0 or NA
		fit_stats <- tryCatch({
			if (sum(!is.na(regTab$y)) < 2 || length(unique(regTab[[target_col]])) < 2) {
				stop("Insufficient data")
			}

			fit <- lm(test_formula, data = regTab)
			sum_fit <- summary(fit)$coefficients

			coeff_idx <- which(grepl(target_col, rownames(sum_fit)))[1]

			if (is.na(coeff_idx)) stop("Term not found")

			sum_fit[coeff_idx, c(1, 4)] # Estimate (log2FC-ish) and p-value

		}, error = function(e) {
			return(c(Estimate = NA, `Pr(>|t|)` = NA))
		})

		return(fit_stats)
	})


	res <- as.data.frame(do.call(rbind, resList))
	colnames(res) <- c("log2FoldChange", "pval")
	rownames(res) <- rownames(exprData)
	res$baseMean <- rowMeans(sub_expr, na.rm = TRUE)
	res$padj <- p.adjust(res$pval, method = "BH")
	return(res[, c("baseMean", "log2FoldChange", "pval", "padj")])
}

checkGSDSobj <- function(x, varname){
	if(!inherits(x, "SummarizedExperiment")) stop(varname, " must be a SummarizedExperiment produced by `computeActivationScore` or `GSDS`")
	if(0 %in% dim(rowData(x))) stop(varname, " must be a SummarizedExperiment with rowData produced by `computeActivationScore` or `GSDS`")
	if(!"activScoreMat" %in% names(assays(x))) stop(varname, " must be a SummarizedExperiment with an `activScoreMat` assay produced by `computeActivationScore` or `GSDS`")
	return(TRUE)
}


checkDB <- function(db) {
	if (!is.list(db)) {
		stop("DBsets must be a list")
	}
	all_elements_are_chars <- vapply(db, function(sub_db) {
		if (!is.list(sub_db)) {
			return(FALSE)
		}
		all_chars <- vapply(sub_db, function(x) {
			is.character(x)
		}, logical(1))

		all(all_chars)
	}, logical(1))

	# If all elements are characters, the sum should be equal to the length of the db.
	all(all_elements_are_chars)
}


getORAstats <- function(intersectionSize, setSizes, universeSize){
	nSet <- length(setSizes)
	expected <- prod(setSizes) / (universeSize ^ (nSet - 1))
	log2OR = NULL
	if (nSet == 2) {
		pval <-
			phyper(
				q = intersectionSize - 0.5,
				m = setSizes[1],
				n = universeSize - setSizes[1],
				k = setSizes[2],
				lower.tail = FALSE
			)
		n10 <- setSizes[1] - intersectionSize
		n01 <- setSizes[2] - intersectionSize
		n00 <- universeSize - n10 - n01 - intersectionSize
		log2OR <- log2( (intersectionSize * n00) / (n10 * n01))
	} else if (nSet > 2) {

		pval <-
			pbinom(intersectionSize - 1,
						 universeSize,
						 expected / universeSize,
						 lower.tail = FALSE)
	} else{
		stop("Number of set should be superior to 1")
	}
	OEdeviation <- (intersectionSize - expected) / sqrt(universeSize)
	res <- list(
		"observed" = intersectionSize,
		"expected" = expected,
		"OEdeviation" = OEdeviation,
		"log2OR" = log2OR,
		"pval" = pval
	)
	if(!is.null(log2OR)) res$log2OR <-log2OR
	return(res)
}

getGeneMarkerStat<-function(gene, designMat, n1, n2){
	aurocRes <- aurocCPP(
		score = gene,
		boolVect = designMat[, 2],
		n1 = n1,
		n2 = n2
	)
	lmOut <- .lm.fit(designMat, gene)
	pval <- pvalLmFit(
		lmOut$residuals,
		lmOut$coefficients,
		p = lmOut$rank,
		qr = lmOut$qr
	)[2]
	coef <- lmOut$coefficients[2]
	score <-
		sign(coef) * prod(c(abs(coef), abs(aurocRes - 0.5),
												min(-log10(pval), 324))) ^ (1 / 3)
	#min(-log10(pval),324): avoid Inf
	return(c(coef, aurocRes, score, pval))
}

#' Compute a dataframe with marker metrics describing best gene marker per group
#' of samples.
#'
#' @param BPPARAM A BPPARAM object as return by [BiocParallel::bpparam()]. Used
#'   for multi-threading.
#' @param returnAsList Return a list where each element is dataframe containing
#'   the marker metrics of a group.
#' @param data A dataframe with genes as rows and samples as columns.
#'   Can also be a `SummarizedExperiment` or `SummarizedExperiment` object.
#' @param groups A vector of group names, same size as the number of columns in
#'   `data`.
#' @param transpose If TRUE, the input data is transposed before processing.
#'   Default is TRUE (feature as rows, samples as columns).
#' @param se_assay Integer or character, if `data` is a
#'   `SummarizedExperiment` related object, the assay name to use.
#'
#' @return A dataframe containing four column per group: Log2(Fold-Change),
#'   AUROC, marker score (see details), p-value and BH adjusted p-value.
#'   If `data` is a `SummarizedExperiment` related object and `returnAsList` is
#'   `FALSE`, the function will add the marker metrics to `rowData`.
#' @details LogFC and pvalues are computed from a linear modelling of the data.
#'
#' Score is consisting of the geometrical mean of absolute LogFC, absolute(auroc
#' - 0.5), and -log10(pval), then signed by the logFC: score = sign(logFC) ×
#' gmean( abs(logFC), abs(aurocRes-0.5), -log10(pval) )
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' markerData <- getMarkers(bulkLogCounts,sampleAnnot$culture_media)
#' sce <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = bulkLogCounts))
#' sce <- getMarkers(sce, sampleAnnot$culture_media, se_assay = "counts")
#' rowData(sce) |> head()
getMarkers <- function (data,
												groups,
												transpose = TRUE,
												BPPARAM = BiocParallel::SerialParam(),
												returnAsList = FALSE,
												se_assay = "logcounts") {
	if (is.null(BPPARAM))
		BPPARAM <- BiocParallel::bpparam()

	sce_obj <-NULL
	if (inherits(data, "SummarizedExperiment")) {
		sce_obj <- data
		data <- assay(sce_obj, se_assay)
	}

	if(!transpose) data <- t(data)
	groups <- as.factor(make.names(groups))
	grpLvl <- levels(groups)
	resDF <- BiocParallel::bplapply(grpLvl,BPPARAM = BPPARAM,
																	FUN = function(group) {
																		logicGroup <- groups == group
																		designMat <- cbind(1, logicGroup)
																		n1 <- as.integer(sum(!logicGroup))
																		n2 <- as.integer(sum(logicGroup))
																		resDFgroup <- apply(data, 1,
																												getGeneMarkerStat, designMat, n1, n2) |>
																			t() |> data.frame()
																		colnames(resDFgroup) <-
																			c("lfc", "auroc", "score", "pval")
																		resDFgroup$padj <-
																			p.adjust(resDFgroup$pval, method = "BH")
																		resDFgroup
																	})
	names(resDF) <- grpLvl

	if (returnAsList) return(resDF)
	resDF<-do.call("cbind", resDF)
	if(is.null(sce_obj)){
		return(resDF)
	}else{
		rowData(sce_obj) <- cbind(rowData(sce_obj),resDF)
		return(sce_obj)
	}

}



#' Extract a specific feature/metric (pval, logFC...) from a marker result
#' dataframe
#'
#' @param markerData A dataframe returned by [oob::getMarkers].
#'   Can also be a `SummarizedExperiment` or `SummarizedExperiment` object where
#'   were [oob::getMarkers] has been performed.
#' @param feature The name of the feature that has to be extracted.
#'
#' @return
#' A matrix containing only the wanted feature where each column is a group.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' markerData <- getMarkers(bulkLogCounts,sampleAnnot$culture_media,
#'     BPPARAM=BiocParallel::SerialParam())
#' extractFeatureMarkerData(markerData) |> head()
#' sc <- SummarizedExperiment(assays = list(counts = bulkLogCounts),
#'    colData = sampleAnnot,
#'    rowData = markerData)
#' extractFeatureMarkerData(sc) |> head()
extractFeatureMarkerData <-
	function(markerData, feature = "score") {
		if (inherits(markerData, "SummarizedExperiment")) {
			markerData <- rowData(markerData)
		}
		featureStrLen <- nchar(feature) + 1
		columns2Pick <-
			grep(paste0("^.*\\.", feature, "$"),
					 colnames(markerData),
					 value = TRUE)
		markerData <- markerData[, columns2Pick, drop = FALSE] |> as.matrix()
		grpName <- colnames(markerData)
		colnames(markerData) <-
			substr(grpName, 1, nchar(grpName) - featureStrLen)
		return(markerData)
	}



#' Complex Heatmap wrapper optimized for RNA-Seq analyses...
#'
#' @inheritParams ComplexHeatmap::Heatmap
#' @param matrix A matrix. Either numeric or character.
#'   If it is a simple vector, it will be converted to a one-column matrix.
#'   Can also be a `SummarizedExperiment` or `SummarizedExperiment` object.
#' @param preset A value from `"expr"`, `"cor"`, `"dist"` or `NULL`. Change
#'   other arguments given a specific preset (default preset if NULL).
#' @param row_name_autofontsize Logical, should row names font size automatically
#'   adjusted to the number of row?
#' @param column_name_autofontsize Logical, should column names font size
#'   automatically adjusted to the number of columns?
#' @param scale Logical. Divide rows of `matrix` by their standard deviation. If
#'   NULL determined by preset.
#' @param center Logical. Subtract rows of `matrix` by their average. If NULL
#'   determined by preset.
#' @param return_heatmap Logical, return the plot as a Heatmap object or print it
#'   in the current graphical device.
#' @param row_names_additionalgpar List. Additional parameter passed to `gpar`
#'   for row names.
#' @param column_names_additionalgpar List. Additional parameter passed to `gpar`
#'   for column names.
#' @param row_names_fontsize Numeric. Font size of the row names
#' @param column_names_fontsize Numeric. Font size of the column names
#' @param border Logical. Whether draw border. The value can be logical or a
#'   string of color.
#' @param colorScale A vector of colors that will be used for mapping colors to
#'   the main heatmap.
#' @param colorScaleFun A function that map values to colors. Used for the main
#'   heatmap. If not NULL this will supersede the use of the `colorScale`
#'   argument.
#' @param midColorIs0 Logical. Force that 0 is the midColor.  If NULL turned on
#'   if the matr.
#' @param probs A numeric vector (between 0 and 1) same length as color or NULL.
#'   Quantile probability of the values that will be mapped to colors.
#' @param useProb Logical. Use quantile probability to map the colors. Else the
#'   min and max of values will be mapped to first and last color and
#'   interpolated continuously.
#' @param minProb A numeric value (between 0 and 1). If `useProb=TRUE` and
#'   `probs=NULL` this will be the quantile of the value for the first color,
#'   quantile will be mapped continuously as to the maxProb.
#' @param maxProb A numeric value (between 0 and 1).
#' @param colData A vector of factor, character, numeric or logical. Or, a
#'   dataframe of any of these type of value. The annotation that will be
#'   displayed on the heatmap.
#' @param colorAnnot List or NULL. Precomputed color scales for the `colData`.
#'   Color scales will be only generated for the features not described. Must be
#'   in the format of a list named by columns of `annots`. Each element contains
#'   the colors at breaks for continuous values. In the case of factors, the
#'   colors are named to their corresponding level or in the order of the
#'   levels.
#' @param show_grid Logical. Draw a border of each individual square on the
#'   heatmap. If NULL automatically true if number of values < 500.
#' @param gpar_grid Gpar object of the heatmap grid if `show_grid`.
#' @param show_values Logical. Show values from the matrix in the middle of each
#'   square of the heatmap.
#' @param nsignif Integer. Number of significant digits showed if `show_values`.
#' @param is_square_heatmap Logical or NULL. Apply clustering columns on rows. If NULL
#'   automatically turned TRUE if `ncol==nrow` and col/rownames are the same.
#' @param se_assay Integer or character, if `data`
#'   is a `SummarizedExperiment` related object, the assay name to use.
#' @param ... Other parameters passed to `Heatmap`.
#'
#' @return A Heatmap object if `return_heatmap` or print the Heatmap in the
#'   current graphical device.
#' @export
#'
#' @seealso [genTopAnnot()], [genRowAnnot()]
#'
#' @details
#'
#' A preset attributes a list of default values for each argument. However, even
#' if a preset is selected, arguments precised by the user precede the preset. #
#' Default arguments ## preset is `NULL`
#' ```
#' clustering_distance_rows = covDist #see covDist for more details
#' clustering_distance_columns = covDist
#' name="matrix"
#' colorScale=c("#2E3672","#4B9AD5","white","#FAB517","#E5261D")
#' center=TRUE
#' scale=FALSE
#' ```
#' ## preset is `"expr"` (expression)
#' ```
#' clustering_distance_rows = covDist
#' clustering_distance_columns = covDist
#' name="centered log expression"
#' colorScale=    c("darkblue","white","red2")
#' row_names_additionalgpar=list(fontface="italic")
#' center=TRUE
#' scale=FALSE
#' ```
#' ## preset is `"cor"` (correlation)
#' ```
#' clustering_distance_rows ="euclidean"
#' clustering_distance_columns ="euclidean"
#' name="Pearson correlation"
#' colorScale=c("darkblue","white","#FFAA00")
#' center=FALSE
#' scale=FALSE
#' ```
#' ## preset is `"dist"` (distance)
#' ```
#' clustering_distance_rows ="euclidean"
#' name="Euclidean distance"
#' colorScale=c("white","yellow","red","purple")
#' center=FALSE
#' scale=FALSE
#' ```
#'
#' ## preset is `"vanilla"` (don't transform value, same as default
#' ComplexHeatmap)
#' ```
#' clustering_distance_rows ="euclidean"
#' name="matrix"
#' colorScale=c("#2E3672","#4B9AD5","white","#FAB517","#E5261D")
#' center=FALSE
#' scale=FALSE
#' ```
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("DEgenesPrime_Naive")
#'
#' bestDE <- rownames(DEgenesPrime_Naive)[whichTop(DEgenesPrime_Naive$pvalue,
#'                                           decreasing = FALSE,
#'                                           top = 50)]
#' HEATmap(
#'     matrix(rnorm(50), ncol = 5),
#'     preset = NULL,
#'     show_values = TRUE,
#'     nsignif = 2
#' )
#'
#' HEATmap(bulkLogCounts[bestDE, ],
#'   colData = sampleAnnot[, c("culture_media", "line")])
#' HEATmap(
#' 	bulkLogCounts[bestDE[seq_len(5)],
#' 	   rownames(sampleAnnot)[sampleAnnot$culture_media %in%
#'     c("T2iLGO", "KSR+FGF2")]]
#' )
#'
#' corDat <- cor(bulkLogCounts)
#' HEATmap(corDat, preset = "cor")
#' HEATmap(
#'     corDat,
#'     preset = "cor",
#'     center = TRUE,
#'     colorScaleFun = circlize::colorRamp2(c(-0.2, 0, 0.2),
#'       c("blue", "white", "red"))
#' )
#' se <- SummarizedExperiment(assays = list(logcounts = bulkLogCounts),
#'    colData = sampleAnnot)
#' HEATmap(se[bestDE[seq_len(5)],], colData = c("line", "culture_media"))
HEATmap <-
	function(matrix,
					 preset = c("expr", "cor", "dist", "vanilla",NULL),
					 clustering_distance_rows = NULL,
					 clustering_distance_columns = NULL,
					 clustering_method_columns = "ward.D2",
					 clustering_method_rows = "ward.D2",
					 row_name_autofontsize = TRUE,
					 column_name_autofontsize = TRUE,
					 row_names_fontsize = NULL,
					 column_names_fontsize = NULL,
					 scale = NULL,
					 center = NULL,
					 return_heatmap = FALSE,
					 name = NULL,
					 row_names_additionalgpar = list(),
					 column_names_additionalgpar = list(),
					 border = TRUE,
					 colorScale = NULL,
					 colorScaleFun = NULL,
					 midColorIs0= NULL,
					 probs = NULL,
					 useProb = TRUE,
					 minProb = 0.05,
					 maxProb = 0.95,
					 cluster_rows = NULL,
					 cluster_columns = NULL,
					 colData = NULL,
					 colorAnnot = NULL,
					 show_grid = NULL,
					 gpar_grid = gpar(col = "black"),
					 show_values = FALSE,
					 nsignif = 3,
					 column_dend_reorder = FALSE,
					 row_dend_reorder = FALSE,
					 is_square_heatmap = NULL,
					 row_split = NULL,
					 column_split = NULL,
					 se_assay = "logcounts",
					 ...) {

		if (inherits(matrix, "SummarizedExperiment")) {
			if (!is.null(colData)) {
				if(inherits(colData, "character")){
					if (sum(!colData %in% colnames(colData(matrix))==0)) {
						colData <- data.frame(colData(matrix)[colData])
					} else {
						stop(
							setdiff(colData, colnames(colData(matrix))),
							"not found in colData of SummarizedExperiment"
						)
					}
				}
			}
			matrix <- assay(matrix, se_assay)
		}
		args <- list()
		preset <- preset[1]
		if (is.null(preset)) {
			if (is.null(clustering_distance_rows))
				clustering_distance_rows <- cosineDist
			if (is.null(clustering_distance_columns))
				clustering_distance_columns <- cosineDist
			if (is.null(name))
				name <- "matrix"
			if (is.null(colorScale))
				colorScale <- c("#2E3672", "#4B9AD5", "white", "#FAB517", "#E5261D")
			if (is.null(row_names_additionalgpar))
				row_names_additionalgpar <- list()
			if (is.null(center))
				center <- TRUE
			if (is.null(scale))
				scale <- FALSE
		} else if (preset == "expr") {
			if (is.null(clustering_distance_rows))
				clustering_distance_rows <- cosineDist
			if (is.null(clustering_distance_columns))
				clustering_distance_columns <- cosineDist
			if (is.null(name))
				name <- "centered log expression"
			if (is.null(colorScale))
				colorScale <- c("darkblue", "white", "red2")
			if (is.null(row_names_additionalgpar))
				row_names_additionalgpar <- list(fontface = "italic")
			if (is.null(center))
				center <- TRUE
			if (is.null(scale))
				scale <- FALSE
		} else if (preset == "cor") {
			if (is.null(clustering_distance_rows))
				clustering_distance_rows <- "euclidean"
			if (is.null(clustering_distance_columns))
				clustering_distance_columns <- "euclidean"
			if (is.null(name))
				name <- "Pearson\ncorrelation"
			if (is.null(colorScale))
				colorScale <- c("darkblue", "white", "#FFAA00")
			if (is.null(row_names_additionalgpar))
				row_names_additionalgpar <- list()
			if (is.null(center))
				center <- FALSE
			if (is.null(scale))
				scale <- FALSE
		} else if (preset == "dist") {
			if (is.null(clustering_distance_rows))
				clustering_distance_rows <- "euclidean"
			if (is.null(clustering_distance_columns))
				clustering_distance_columns <- "euclidean"
			if (is.null(name))
				name <- "Euclidean\ndistance"
			if (is.null(colorScale))
				colorScale <- c("white", "yellow", "red", "purple")
			if (is.null(row_names_additionalgpar))
				row_names_additionalgpar <- list()
			if (is.null(center))
				center <- FALSE
			if (is.null(scale))
				scale <- FALSE
		} else if (preset == "vanilla") {
			if (is.null(clustering_distance_rows))
				clustering_distance_rows <- cosineDist
			if (is.null(clustering_distance_columns))
				clustering_distance_columns <- cosineDist
			if (is.null(name))
				name <- "matrix"
			if (is.null(colorScale))
				colorScale <- c("#2E3672", "#4B9AD5", "white", "#FAB517", "#E5261D")
			if (is.null(row_names_additionalgpar))
				row_names_additionalgpar <- list()
			if (is.null(center))
				center <- FALSE
			if (is.null(scale))
				scale <- FALSE
		} else{
			stop("preset must equal to one of this value: NULL, ",
					 "'expr', 'cor', 'dist', 'vanilla'")
		}
		matrix <- as.matrix(matrix)

		if (min(apply(matrix, 1, sd, na.rm = TRUE)) == 0 &
				(scale |
				 identical(corrDist, clustering_distance_rows))) {
			warning(
				"some row have a 0 sd. sd-based method ",
				"(correlation distance, scaling) ",
				"will be deactivated or switched."
			)
			scale <- FALSE
			if (identical(corrDist, clustering_distance_rows)) {
				args$clustering_distance_rows <- "euclidean"
			}
		}
		if (scale |
				center)
			matrix <-
			rowScale(matrix, scaled = scale, center = center)
		if (is.null(midColorIs0)) {
			if (min(matrix, na.rm = TRUE) < 0 & max(matrix, na.rm = TRUE) > 0) {
				midColorIs0<- TRUE
			} else{
				midColorIs0<- FALSE
			}
		}
		if (is.null(is_square_heatmap)) {
			if (nrow(matrix) == ncol(matrix) &
					identical(colnames(matrix),rownames(matrix))) {
				is_square_heatmap <- TRUE
				warning("colnames and rownames are identical, ",
								"is_square_heatmap is set to TRUE")
			} else{
				is_square_heatmap <- FALSE
			}
		}
		if (is_square_heatmap) {
			if (is.null(cluster_columns)) {
				cluster_columns <-
					hierarchicalClustering(
						matrix,
						transpose = FALSE,
						method.dist = clustering_distance_columns,
						method.hclust = clustering_method_columns
					)
			}
			args$cluster_rows <- cluster_columns
			args$cluster_columns <- cluster_columns
		} else{
			if (is.null(cluster_rows)) {
				args$clustering_method_rows <- clustering_method_rows
				args$clustering_distance_rows <-
					clustering_distance_rows
			} else{
				args$cluster_rows <- cluster_rows
			}
			if (is.null(cluster_columns)) {
				args$clustering_method_columns <- clustering_method_columns
				args$clustering_distance_columns <-
					clustering_distance_columns
			} else{
				args$cluster_columns <- cluster_columns
			}
		}

		if (is.null(colorScaleFun)) {
			colorScaleFun <-
				computeColorScaleFun(
					colors = colorScale,
					values = unlist(matrix),
					useProb = useProb,
					probs = probs,
					minProb = minProb,
					maxProb = maxProb,
					midColorIs0= midColorIs0,
					returnColorFun = TRUE
				)
		}
		args$col <- colorScaleFun
		if (is.null(show_grid)) {
			if (nrow(matrix) * ncol(matrix) < 500) {
				show_grid <- TRUE
			} else{
				show_grid <- FALSE
			}
		}
		if (show_grid) {
			args$rect_gp <- gpar_grid
		}
		if (show_values) {
			args$cell_fun <- function(j, i, x, y, w, h, col) {
				#dark or light background .
				if (colSums(col2rgb(col)) < 382.5)
					col <- "white"
				else
					col <- "black"
				grid.text(
					as.character(
						signif(matrix[i, j], nsignif)
					), x, y, gp = gpar(col = col)
				)
			}
		}
		args$row_names_gp <- row_names_additionalgpar
		args$column_names_gp <- column_names_additionalgpar
		if(!is.null(row_names_fontsize)){
			args$row_names_gp$fontsize <- row_names_fontsize
		} else if (row_name_autofontsize){
			args$row_names_gp <- do.call("autoGparFontSizeMatrix",
																	 c(list(nrow(matrix)), row_names_additionalgpar))
		}
		if(!is.null(column_names_fontsize)){
			args$column_names_gp$fontsize <- column_names_fontsize
		} else if (column_name_autofontsize){
			args$column_names_gp <- do.call("autoGparFontSizeMatrix",
					c(list(ncol(matrix)), column_names_additionalgpar))
		}
		if (!is.null(colData)) {
			args$top_annotation <- genTopAnnot(colData, colorAnnot)
		}

		args$column_dend_reorder <-
			column_dend_reorder
		args$row_dend_reorder <- row_dend_reorder
		args$row_split <- row_split
		args$column_split <- column_split
		args$matrix <- matrix
		args$name <- name
		args$border <- border
		args <- c(args, list(...))

		ht <- do.call("Heatmap", args)
		if (return_heatmap) {
			return(ht)
		} else{
			print(ht)
		}
	}




#' Plot best marker per group on a tidy Heatmap
#'
#' @param countTable A matrix of numeric with samples as columns (in the RNA-Seq
#'   context, log counts)
#' @param group A feature of factor/character, same length as number of sample.
#'   Describe group of each sample (for example clusters).
#' @param markerData A matrix describing marker scores for each group.
#' @param topn Number of marker to plot, rankes from the best
#' @param show_column_names Whether show column names.
#' @param ... Arguments passed to `HEATmap`
#'
#' @inheritParams HEATmap
#' @inheritParams drawSamplePerGroup
#'
#' @details Draw the same number of observation from each condition/group, and
#'   take the top n marker per group.
#'   The Heatmap is sliced by group of samples, and by their respective markers.
#'
#' @return A Heatmap object if `return_heatmap` or print the Heatmap in the
#'   current graphical device.
#' @export
#'
#' @seealso [HEATmap()]
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' markerData <- getMarkers(bulkLogCounts,sampleAnnot$culture_media,
#'     BPPARAM=BiocParallel::SerialParam())
#' htMarker(bulkLogCounts,  group=sampleAnnot$culture_media,
#'     markerData=extractFeatureMarkerData(markerData),
#'     colData=sampleAnnot[c("line","passage")])

htMarker <-
	function(countTable,
					 group,
					 markerData,
					 colData = NULL,
					 topn = 5,
					 maxDrawSize = NULL,
					 minDrawSize = NULL,
					 replace = FALSE,
					 return_heatmap = FALSE,
					 show_column_names = FALSE,
					 ...)
	{
		if (length(group) != ncol(countTable))
			stop("Length of group should be the same as number of columns in ",
					 "countTable")
		group <- as.factor(make.names(group))
		if (sum(!sort(colnames(markerData)) == sort(levels(group))) >
				0)
			stop("colnames of markerData must ",
					 "correspond to the levels of group")
		if(is.null(colnames(countTable)))
			colnames(countTable) <- make.names(
				paste0("X",seq_len(ncol(countTable))))
		names(group) <- colnames(countTable)
		if(!all(rownames(markerData) %in% rownames(countTable))) {
			stop("rownames of markerData must ",
					 "exist in the rownames of countTable")
		}

		drawPerGroup<-drawSamplePerGroup(names(group),group,
																		 minDrawSize = maxDrawSize,
																		 maxDrawSize = maxDrawSize, replace = replace)
		drawCells <- names(drawPerGroup)

		topMarker <- apply(markerData, 2, function(x) {
			names(x)[order(x, decreasing = TRUE)][seq_len(topn)]
		})
		topMarker <- as.list(data.frame(topMarker))
		HEATmap(
			countTable[unlist(topMarker), drawCells],
			colData = colData[drawCells,],
			column_split = drawPerGroup,
			row_split = VectorListToFactor(topMarker),
			cluster_row_slices = FALSE,
			cluster_column_slices = FALSE,
			return_heatmap = return_heatmap,
			show_column_names = show_column_names,
			...
		)
	}


#' Generate a top annotation for ComplexHeatmap
#'
#' @param annot A vector of factor, character, numeric or logical. Or, a
#'   dataframe of any of these type of value. The annotation that will be
#'   displayed on the heatmap.
#' @param colorScales List or NULL. Precomputed color scales. Color scales will
#'   be only generated for the features not described. Must be in the format of
#'   a list named by columns of `annots`. Each element contains the colors at
#'   breaks for continuous values. In the case of factors, the colors are named
#'   to their corresponding level or in the order of the levels.
#' @param border Logical. Whether draw border. The value can be logical or a
#'   string of color.
#' @param ... Other parameters passed to `genColorsForAnnots`.
#'
#' @return A HeatmapAnnotation object. Can be used for example in the
#'   `top_annotation` argument of `Heatmap`.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("DEgenesPrime_Naive")
#'
#' bestDE <-
#'     rownames(DEgenesPrime_Naive)[
#'         whichTop(DEgenesPrime_Naive$pvalue,
#'         decreasing = FALSE,
#'         top = 50)
#'     ]
#' ComplexHeatmap::Heatmap(rowScale(bulkLogCounts[bestDE, ]), top_annotation =
#'     genTopAnnot(sampleAnnot[, c("culture_media", "line")]))
genTopAnnot <- function(annot,
												colorScales = NULL,
												border = TRUE,
												...) {
	if (is.factor(annot) | is.data.frame(annot))
		annot <- droplevels(annot)
	if (is.vector(annot) | is.factor(annot)) {
		annot <- data.frame(Annotation = annot)
		if (is.list(colorScales))
			colnames(annot) <- names(colorScales)[1]
		if ((!is.null(colorScales)) &
				!is.list(colorScales))
			colorScales <- list("Annotation" = colorScales)
	}
	colorScales <-
		genColorsForAnnots(
			annots = annot,
			colorScales = colorScales,
			returnContinuousFun = TRUE,
			...
		)
	HeatmapAnnotation(df = annot,
										col = colorScales,
										border = border)
}


#' Generate a row annotation for ComplexHeatmap
#'
#' @param annot A vector of factor, character, numeric or logical. Or, a
#'   dataframe of any of these type of value. The annotation that will be
#'   displayed on the heatmap.
#' @param colorScales List or NULL. Precomputed color scales. Color scales will
#'   be only generated for the features not described. Must be in the format of
#'   a list named by columns of `annots`. Each element contains the colors at
#'   breaks for continuous values. In the case of factors, the colors are named
#'   to their corresponding level or in the order of the levels.
#' @param border Logical. Whether draw border. The value can be logical or a
#'   string of color.
#' @param ... Other parameters passed to `genColorsForAnnots`.
#'
#' @return A HeatmapAnnotation object. Can be used for example in the
#'   `top_annotation` argument of `Heatmap`.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("DEgenesPrime_Naive")
#'
#' library(ComplexHeatmap)
#'
#' bestDE <- rownames(DEgenesPrime_Naive)[whichTop(DEgenesPrime_Naive$pvalue,
#'                                           decreasing = FALSE,
#'                                           top = 50)]
#' Heatmap(rowScale(bulkLogCounts[bestDE, ]) |> t(),
#'         right_annotation  = genRowAnnot(
#'           sampleAnnot[, c("culture_media", "line")])
#'         )
genRowAnnot <-
	function(annot,
					 colorScales = NULL,
					 border = TRUE,
					 ...) {
		if (is.factor(annot) | is.data.frame(annot))
			annot <- droplevels(annot)
		if (is.vector(annot) | is.factor(annot)) {
			annot <- data.frame(Annotation = annot)
			if (is.list(colorScales))
				colnames(annot) <- names(colorScales)[1]
			if ((!is.null(colorScales)) & !is.list(colorScales))
				colorScales <- list(Annotation = colorScales)
		}
		colorScales <-
			genColorsForAnnots(
				annots = annot,
				colorScales = colorScales,
				returnContinuousFun = TRUE,
				...
			)
		ComplexHeatmap::rowAnnotation(df = annot,
																	col = colorScales,
																	border = border)
	}




#' Displaying most neg and pos gene contribution on a GSDS heatmap.
#'
#' @param contributions A list named by gene set Contains vector of gene contribution to activation scores, named by genes.
#' @param maxgene_shownperside Integer. Number of best pos/neg rank displayed on the annotation
#' @param fontface Character. Fontface of the genes text
#' @param fontsize Unit. Fontsize of the genes text. If NULL computed automatically
#' @param width Unit. Width of the annotation heatmap. If NULL computed automatically
#' @param max_width Unit. If width is computed automatically, maximum width of the contribution annotation
#' @param additional_gp Additional gpar argument passed as a list, see ?grid::gpar
#'
#' @return A HeatmapAnnotation object. Genes on the left side of the vertical line are contributing negatively to the activation score, and positively on the right side.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("keggHuman")
#'
#' resGSDS <- GSDS(
#' 	data = bulkLogCounts,
#' 	DBsets = keggHuman,
#' 	colData = sampleAnnot,
#' 	contrast = c("culture_media", "T2iLGO", "KSR+FGF2"),
#' 	return_seobject = TRUE
#' )
#' resGSDS.sel <- resGSDS[abs(rowData(resGSDS)$log2FoldChange) > 10]
#'
#' nrow(resGSDS.sel)
#'
#' HEATmap(
#' 		assay(resGSDS.sel),
#'     midColorIs0= TRUE,
#'     center = FALSE,
#'     name = "Activation score",
#'     preset = NULL,
#'     colData = sampleAnnot["culture_media"],
#'     right_annotation = ComplexHeatmap::rowAnnotation(
#'         "gene contribution" =
#'             GSDS.HeatmapAnnot(
#'                 contributions = rowData(resGSDS)$contributions,
#'                 fontsize = unit(8,"points")
#'
#'             )
#'     ),
#'     row_names_side = "left",
#'     row_dend_side = "left",
#'     row_names_max_width = unit(340, "points")
#' )
GSDS.HeatmapAnnot <-
	function(contributions,
					 maxgene_shownperside = 3,
					 width = NULL,
					 max_width = unit(500, "points"),
					 fontface = "italic",
					 fontsize = NULL,
					 additional_gp=list()) {

		if(!is.list(contributions)) stop("contributions must be a list")
		n_gs <- length(contributions)
		if(is.null(fontsize)) fontsize <- unit(autoGparFontSizeMatrix(n_gs)$fontsize, "points")
		if(is.null(width)) width <- max(unit(50, "points"), min(35 * fontsize, max_width))
		if(!is.unit(width)) stop("width must be a unity create by grid::unit")
		if(!is.unit(fontsize)) stop("fontsize must be a unity create by grid::unit")
		ComplexHeatmap::AnnotationFunction(
			fun = function(index, k, n) {
				generateGeneNameSide <- function(text,x, y, just,fontsize, fontface = "italic",additional_gp = NULL){
					grid.text(
						paste0(text, collapse = "  "),
						just = just,
						x = x,
						y,
						default.units = "native",
						gp = do.call("gpar",c(list("fontsize" = fontsize, "fontface" = fontface),additional_gp))
					)
				}
				pushViewport(viewport(
					xscale = c(0, 10),
					yscale = c(0.5, length(index) + 0.5)
				))
				grid.lines(.5, c(0, 1))
				i <- length(index)
				for (sel in index) {
					contrib <- sort(contributions[[sel]])

					left <- names(contrib)[contrib < 0]
					right <- names(contrib)[contrib > 0]
					left <- left[1:min(maxgene_shownperside, length(left))]
					right <-
						right[max(1, (length(right) - maxgene_shownperside) + 1):length(right)]
					if (!is.null(left[1]) & !is.na(left[1])) {
						generateGeneNameSide(
							left,
							x = 4.5,
							y = i,just = "right",
							fontsize = fontsize,
							fontface = fontface,
							additional_gp = additional_gp
						)
					}
					if (!is.null(right[1]) & !is.na(right[1])) {
						generateGeneNameSide(
							right,
							x = 5.5,
							y = i,just = "left",
							fontsize = fontsize,
							fontface = fontface,
							additional_gp = additional_gp
						)
					}
					i <- i - 1
				}
				popViewport()
			},
			var_import = list(
				"contributions" = contributions,
				"maxgene_shownperside" = maxgene_shownperside,
				"width" = width,
				"fontface" = fontface,
				"fontsize" = fontsize,
				"additional_gp" = additional_gp
			),
			subsettable = FALSE,
			width = width,
			which = "row"
		)
	}




#' Plot the activation scores and the top gene contributions in one heatmap
#' @inheritParams HEATmap
#' @param GSDS.object A `SummarizedExperiment` as returned by `computeActivationScore` or `GSDS`
#' @param row_names_max_width Unit. Max width of the heatmap gene set
#' @param contrib_max_width Unit. If width is computed automatically, maximum width of the contribution annotation
#' @param maxgene_shownperside Integer. Number of best pos/neg rank displayed on the annotation
#' @param contrib_width Unit. Width of the gene contribution annotation. If NULL computed automatically
#' @param contrib_fontface Character. Fontface of the genes contribution annotation
#' @param contrib_fontsize  Unit. Fontsize of the genes contribution annotation. If NULL computed automatically
#' @param contrib_additional_gp Additional gpar  argument passed as a list to the contribution annotation, see ?grid::gpar
#' @param ... Arguments passed to HEATmap
#'
#' @returns A heatmap of the activation scores, names by gene sets, and with the top gene contributions on the right.
#' a selection of genes with a positive or negative contribution to the activation score is represented each set.
#' Positive contribution (at the right of the vertical bar), means that up regulation increases the activation score, whereas
#' a negative contribution (left of the bar) means that an increase of expression will decrease the activation score.
#' This makes it easy to interpret whether a positive activation score means that a pathway is "biologically" on or off,
#' or assess the redundance of gene sets.
#' @seealso [HEATmap()], [GSDS.HeatmapAnnot()]
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("keggHuman")
#'
#' resGSDS <- GSDS(
#' 	data = bulkLogCounts,
#' 	DBsets = keggHuman,
#' 	colData = sampleAnnot,
#' 	contrast = c("culture_media", "T2iLGO", "KSR+FGF2"),
#' 	return_seobject = TRUE
#' )
#' resGSDS.sel <- resGSDS[abs(rowData(resGSDS)$log2FoldChange) > 7]
#'
#' GSDS.Heatmap(resGSDS.sel,  colData = sampleAnnot["culture_media"])
#'
GSDS.Heatmap <- function(GSDS.object, colData = NULL, row_names_max_width = unit(340, "points"),
												 maxgene_shownperside = 3,
												 contrib_width = NULL,
												 contrib_max_width = unit(500, "points"),
												 contrib_fontface = "italic",
												 contrib_fontsize = NULL,
												 contrib_additional_gp=list(),...){
	invisible(checkGSDSobj(GSDS.object,"GSDS.object"))
	HEATmap(
		assay(GSDS.object, "activScoreMat"),
		midColorIs0= TRUE,
		center = FALSE,
		name = "Activation score",
		preset = NULL,
		colData =colData,
		right_annotation = rowAnnotation(
			"gene contribution" =
				GSDS.HeatmapAnnot(
					contributions = rowData(GSDS.object)$contributions,
					maxgene_shownperside = maxgene_shownperside,
					width = contrib_width,
					max_width = contrib_max_width,
					fontface = contrib_fontface,
					fontsize = contrib_fontsize,
					additional_gp=contrib_additional_gp
				)
		),
		row_names_side = "left",
		row_dend_side = "left",
		row_names_max_width = row_names_max_width,
		...
	)
}



autoGparFontSizeMatrix <- function(n, ...) {
	n <- max(n, 50)
	n <- min(n, 1000)
	return(gpar(fontsize = 1 / n * 600, ...))
}

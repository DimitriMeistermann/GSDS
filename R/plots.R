

#' Displaying most neg and pos gene contribution on a GSDS heatmap.
#'
#' @param contributions A list named by pathway. Contains vector of gene contribution to activation scores, named by genes.
#' @param maxGeneContribAtOneSide Integer. Number of best pos/neg rank displayed on the annotation
#' @param width Numeric. Width of the annotation heatmap.
#' @param fontsizeFactor Numeric. Font-size of gene names.
#'
#' @return A HeatmapAnnotation object. Genes on the left side of the vertical line are contributing negatively to the activation score, and positively on the right side.
#' @export
#'
#' @examples
#' library(oob, quietly = TRUE)
#'
#' data("bulkLogCounts", package = "oob")
#' data("sampleAnnot", package = "oob")
#'
#' keggDB <- getDBterms(rownames(bulkLogCounts), database = "kegg")
#' geneSetActivScore <-
#'     computeActivationScore(bulkLogCounts, db_terms = keggDB)
#' resGSDS <-
#'     GSDS(
#'         geneSetActivScore = geneSetActivScore,
#'         colData = sampleAnnot,
#'         contrast = c("culture_media", "T2iLGO", "KSR+FGF2"),
#'         db_terms =  keggDB
#'     )
#' bestPathayIndex <- whichTop(resGSDS$padj, top = 30, decreasing = FALSE)
#' selectedActivScores <- t(geneSetActivScore$kegg$activScoreMat[, bestPathayIndex])
#' selectedContributions <- geneSetActivScore$kegg$contribution[bestPathayIndex]
#'
#'
#' heatmap.DM(
#'     selectedActivScores,
#'     midColorIs0 = TRUE,
#'     center = FALSE,
#'     name = "Activation score",
#'     preSet = NULL,
#'     colData = sampleAnnot["culture_media"],
#'     right_annotation = rowAnnotation(
#'         "gene contribution" =
#'             GSDS.HeatmapAnnot(
#'                 contributions = selectedContributions,
#'                 width = unit(12, "cm"),
#'                 fontsizeFactor = 300
#'             )
#'     ),
#'     row_names_side = "left",
#'     row_dend_side = "left",
#'     row_names_max_width = unit(8, "inches"),
#'     autoFontSizeRow = FALSE,
#'     row_names_gp = gpar(fontsize = 1 / length(bestPathayIndex) * 300)
#' )

GSDS.HeatmapAnnot <-
    function(contributions,
             maxGeneContribAtOneSide = 3,
             width = unit(3, "cm"),
             fontsizeFactor = 400) {
        ComplexHeatmap::AnnotationFunction(
            fun = function(index, k, n) {
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
                    left <- left[1:min(maxGeneContribAtOneSide, length(left))]
                    right <-
                        right[max(1, (length(right) - maxGeneContribAtOneSide) + 1):length(right)]
                    if (!is.null(left[1])) {
                        grid.text(
                            paste0(left, collapse = "  "),
                            just = "right",
                            x = 4.5,
                            y = i,
                            default.units = "native",
                            gp = gpar(
                                fontsize = 1 / length(index) * fontsizeFactor,
                                fontface = "italic"
                            )
                        )
                    }
                    if (!is.null(right[1])) {
                        grid.text(
                            paste0(right, collapse = "  "),
                            just = "left",
                            x = 5.5,
                            y = i,
                            default.units = "native",
                            gp = gpar(
                                fontsize = 1 / length(index) * fontsizeFactor,
                                fontface = "italic"
                            )
                        )
                    }
                    i <- i - 1
                }
                popViewport()
            },
            var_import = list(
                "contributions" = contributions,
                "maxGeneContribAtOneSide" = maxGeneContribAtOneSide,
                "width" = width,
                "fontsizeFactor" = fontsizeFactor
            ),
            subsettable = FALSE,
            width = width,
            which = "row"
        )
    }

#' General volcano plot
#'
#' @param d A dataframe containing the data needed for the volcano plot. Must
#'   have column names. Can be `NULL` if `effectSizeCol`, `adjPvalCol` and
#'   `labelCol` are vectors.
#' @param effectSizeCol Column name containing the effect size column
#'   (Log2FoldChange for example). Can also be a vector of numeric containing
#'   the effect size values.
#' @param adjPvalCol Column name containing the adjusted pval column. Can also
#'   be a vector of numeric containing the padj values.
#' @param labelCol Column name containing the feature name column (gene name for
#'   example).  Can also be a vector of character containing the labels.
#' @param padjThres Significativity threshold of adjusted p-value for consider a
#'   feature significant.
#' @param minEffectSize Absolute minimum effect size to consider a feature
#'   significant.
#' @param topShownPerSide Number of feature shown at the left and right side of
#'   the volcano plot.
#' @param returnGraph Logical. Return the graph as a ggplot object instead of
#'   printing it.
#' @param neutralVal Value considered as null effect size.
#' @param ... Parameters passed to geom_repel
#'
#' @return Plot in the current graphical device or a ggplot object if
#'   `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' data(sampleAnnot)
#' volcanoPlot(d = DEgenesPrime_Naive,effectSizeCol = "log2FoldChange",
#'     adjPvalCol = "padj", minEffectSize = 1,
#'     labelCol = rownames(DEgenesPrime_Naive))
#'
volcanoPlot <-
	function(d = NULL,
					 effectSizeCol,
					 adjPvalCol,
					 labelCol,
					 padjThres = 0.05,
					 minEffectSize = 0,
					 topShownPerSide = 15,
					 returnGraph = FALSE,
					 neutralVal = 0,
					 ...) {
		if (is.null(d)) {
			if (!(is.numeric(effectSizeCol) &
						is.numeric(adjPvalCol) &
						length(labelCol) > 1))
				stop("If d is null, other parameters must be vector ",
						 "of the same size")
			d <-
				data.frame(effectSize = effectSizeCol,
									 adjPval = adjPvalCol,
									 label = labelCol)
			effectSizeCol <- "effectSize"
			adjPvalCol <- "adjPval"
			labelCol <- "label"

		} else{
			d <- data.frame(d)
			if (is.null(colnames(d)))
				stop("d must have colnames")
			if (is.numeric(effectSizeCol)) {
				d$effectSize <- effectSizeCol
				effectSizeCol <- "effectSize"
			}
			if (is.numeric(adjPvalCol)) {
				d$adjPval <- adjPvalCol
				adjPvalCol <- "adjPval"
			}
			if (length(labelCol) > 1) {
				d$label <- labelCol
				labelCol <- "label"
			}
			if (!effectSizeCol %in% colnames(d))
				stop(effectSizeCol, "is not a colname of d")
			if (!adjPvalCol %in% colnames(d))
				stop(adjPvalCol, "is not a colname of d")
			if (!labelCol %in% colnames(d))
				stop(labelCol, "is not a colname of d")
		}
		d <-
			d[, c(effectSizeCol, adjPvalCol, labelCol)] |> na.omit()
		scale <-
			scales::trans_new(
				name = "invLog10",
				transform = function(x)
					- log10(x),
				inverse = function(x)
					10 ^ (-x),
				domain = c(0, Inf)
			)
		breaks <- ggplotBreak(d[, adjPvalCol], scale)
		xlims <- max(abs(d[, effectSizeCol] - neutralVal))
		xlims <- c(neutralVal - xlims, neutralVal + xlims)
		isNegEffectSize <- d[, effectSizeCol] < neutralVal
		shownLabel <-
			c(
				d[isNegEffectSize, labelCol][whichTop(
					d[isNegEffectSize, adjPvalCol],
					top = topShownPerSide, decreasing = FALSE
				)],
				d[!isNegEffectSize, labelCol][whichTop(
					d[!isNegEffectSize, adjPvalCol],
					top = topShownPerSide, decreasing = FALSE
				)]
			)
		d$significant <-
			d[, adjPvalCol] < padjThres &
			abs(d[, effectSizeCol] - neutralVal) > minEffectSize
		g <-
			ggplot(d,
						 aes(x = .data[[effectSizeCol]], y = .data[[adjPvalCol]],
						 		color = .data$significant)) +
			scale_y_continuous(trans = scale, breaks = breaks) +
			geom_point() + theme_bw() +
			scale_color_manual(values = c("grey75", "black")) +
			xlim(xlims) +
			ggrepel::geom_text_repel(
				data = d[d[, labelCol] %in% shownLabel,],
				aes(x = .data[[effectSizeCol]], y = .data[[adjPvalCol]],
						label = .data[[labelCol]]),
				inherit.aes = FALSE,
				color = "grey50",
				...
			)
		if (returnGraph) {
			return(g)
		} else{
			print(g)
		}
	}


ggplotBreak<-function(x, scale, m = 5) {
	transValues <- scale$transform(x)
	breaks <-
		labeling::extended(min(transValues), max(transValues), m = m)
	scale$inverse(breaks)
}



#' Compute a color scale function from numeric values by interpolation
#'
#' @param colors A character vector containing the colors.
#' @param values A numeric vector of the value to has to be mapped to colors.
#' @param useProb Logical. Use quantile probability to map the colors. Else the
#'   min and max of values will be mapped to first and last color and
#'   interpolated continuously.
#' @param probs A numeric vector (between 0 and 1) same length as color or NULL.
#'   Quantile probability of the values that will be mapped to colors.
#' @param minProb A numeric value (between 0 and 1). If `useProb=TRUE` and
#'   `probs=NULL` this will be the quantile of the value for the first color,
#'   quantile will be mapped continuously as to the maxProb.
#' @param maxProb A numeric value (between 0 and 1).
#' @param midColorIs0 Logical. Force that 0 return the midColor.
#' @param returnColorFun Logical.Return converted values to colors or the scale
#'   function.
#' @param returnGGscale Logical. Return a ggplot2 gradiantn scale.
#' @param geomAes "fill" or "color". Ggplot layer that will receive the scale.
#' @param geomArgument list of additional argument to pass to the ggplot2
#'   gradiantn scale.
#'
#' @return A vector of colors, or a function if `returnColorFun=TRUE` or a
#'   ggplot scale if `returnGGscale=TRUE`.
#' @export
#'
#' @examples
#' library(ggplot2)
#' values=sort(rnorm(100))
#'
#' scales::show_col(computeColorScaleFun(
#'     colors = c("black", "red"),
#'     values = values,
#'     returnColorFun = FALSE
#' ))
#' scales::show_col(computeColorScaleFun(
#'     colors = c("blue", "white", "red"),
#'     values = values,
#'     returnColorFun = FALSE
#' ))
#' scales::show_col(
#'     computeColorScaleFun(
#'         colors = c("blue", "white", "red"),
#'         values = values,
#'         returnColorFun = FALSE,
#'         midColorIs0 = TRUE
#'     )
#' )
#' scales::show_col(
#'     computeColorScaleFun(
#'         colors = c("blue", "white", "red"),
#'         values = values,
#'         returnColorFun = FALSE,
#'         useProb = TRUE
#'     )
#' )
#' scales::show_col(
#'     computeColorScaleFun(
#'         colors = c("blue", "white", "red"),
#'         values = values,
#'         returnColorFun = FALSE,
#'         useProb = TRUE,
#'         probs = c(.25, .5, .75)
#'     )
#' )
#'
#' colorFun <-
#'     computeColorScaleFun(
#'         colors = c("blue", "white", "red"),
#'         values = values,
#'         returnColorFun = TRUE,
#'         useProb = TRUE
#'     )
#' scales::show_col(c(colorFun(-1), colorFun(0), colorFun(1)))
#'
#' dat <- data.frame(x = rnorm(10),
#'                   y = rnorm(10),
#'                   expr = rnorm(10))
#' ggplot(dat, aes(x = x, y = y, fill = expr)) +
#'     geom_point(size = 5, shape = 21) + theme_bw() +
#'     computeColorScaleFun(
#'         colors = c("blue", "white", "red"),
#'         values = dat$expr,
#'         returnGGscale = TRUE,
#'         useProb = TRUE,
#'         geomAes = "fill"
#'     )
computeColorScaleFun <- function (colors,
																	values,
																	useProb = FALSE,
																	probs = NULL,
																	minProb = 0.05,
																	maxProb = 0.95,
																	midColorIs0 = FALSE,
																	returnColorFun = TRUE,
																	returnGGscale = FALSE,
																	geomAes = "fill",
																	geomArgument = list())
{
	if (is.null(values))
		stop("values cannot be NULL")
	if (!useProb) {
		breaks <- seq(min(values, na.rm = TRUE),
									max(values, na.rm = TRUE),
									length.out = length(colors))
	}
	else {
		if (is.null(probs)) {
			probs <- seq(minProb, maxProb, length.out = length(colors))
		}
		breaks <- quantile(values, probs = probs, na.rm = TRUE)
	}
	if (midColorIs0 & (length(colors) %% 2 == 1)) {
		breaks[ceiling(length(breaks) / 2)] <- 0
	}
	colorFun <- circlize::colorRamp2(breaks = breaks, colors = colors)
	if (returnGGscale) {
		scaledBreaks <- linearScale(values, c(0, 1), returnFunction = TRUE)(breaks)
		if (scaledBreaks[1] > 0) {
			scaledBreaks <- c(0, scaledBreaks)
			colors <- c(colors[1], colors)
		}
		if (scaledBreaks[length(scaledBreaks)] < 1) {
			scaledBreaks <- c(scaledBreaks, 1)
			colors <- c(colors, colors[length(colors)])
		}
		geomArgument$values <- scaledBreaks
		geomArgument$colors <- colors
		return(do.call(paste0("scale_", geomAes, "_gradientn"), geomArgument))
	}
	if (returnColorFun) {
		return(colorFun)
	}
	else {
		return(colorFun(values))
	}
}


#' Generate a list of value/color mapping
#'
#' @param annots Dataframe. Can contain factors or numeric. Contains the values
#'   that has to be mapped.
#' @param colorScales List or NULL. Precomputed color scales. Color scales will
#'   be only generated for the features not described. Must be in the format of
#'   a list named by columns of `annots`. Each element contains the colors at
#'   breaks for continuous values or a mapping function if
#'   `returnContinuousFun=TRUE` (a function that return a color for a given
#'   numeric value). In the case of factors, the colors are named to their
#'   corresponding level, or in the order of the levels.
#' @param discreteFuns A list functions that take a single integer n and return
#'   n colors. If several functions are provided it will be used for each factor
#'   column successively.
#' @param returnContinuousFun Logical. Return a mapping function for continuous
#'   values instead of a vector of colors.
#' @param continuousPalettes A list of color vector. If several vector are
#'   provided it will be used for each numerical column successively.
#' @param ... Parameters passed to `computeColorScaleFun`.
#'
#' @return A list describing the color scale of each column of `annots`, in the
#'   same format than the argument `colorScales`
#' @export
#'
#' @examples
#' data("iris")
#'
#' genColorsForAnnots(iris)
#'
#' precomputedColorScale <-
#'     list(Species = c(
#'         "setosa" = "red",
#'         "versicolor" = "blue",
#'         "virginica" = "grey"
#'     ))
#'
#' genColorsForAnnots(iris, colorScales = precomputedColorScale)
#'
#' colorScales <- genColorsForAnnots(iris, returnContinuousFun = TRUE)
#' colorScales$Sepal.Length(4.5)
#' colorScales$Species
#'
#' library(ComplexHeatmap)
#' Heatmap(
#'     rowScale(t(
#'         iris[, c("Sepal.Length", "Sepal.Width",
#'             "Petal.Length", "Petal.Width")]
#'     ),
#'     center = TRUE, scaled = TRUE),
#'     top_annotation = genTopAnnot(iris["Species"], colorScales =
#'                                      colorScales["Species"])
#' )

genColorsForAnnots <-
	function(annots,
					 colorScales = NULL,
					 discreteFuns = list(gsds_colors, mostDistantColor, mostDistantColor),
					 returnContinuousFun = FALSE ,
					 continuousPalettes = list(
					 	c("#440154", "#6BA75B", "#FDE725"),
					 	c("#2EB538", "#1D1D1B", "#DC0900"),
					 	c("#FFFFC6", "#FF821B", "#950961")
					 ),
					 ...) {
		if (is.null(colnames(annots)))
			stop("annots must have colnames.")
		annotNames <- colorScalesToGen <- colnames(annots)
		newColorScales <- list()
		if (!is.null(colorScales)) {
			for (colorScaleName in names(colorScales)) {
				if (!colorScaleName %in% colorScalesToGen)
					stop(
						"Condition '",
						colorScaleName,
						"' does not match with existing condition names"
					)
				colorScale <- colorScales[[colorScaleName]]
				annotVect <- annots[, colorScaleName]
				if (!is.null(names(colorScale))) {
					#factors
					if (is.numeric(annotVect)) {
						warning(
							colorScaleName,
							" is numeric but encoded as factors (color vector ",
							"has names). It will be converted to factors."
						)
						annots[, colorScaleName] <-
							as.factor(as.character(colData[, colorScaleName]))
						annotVect <- annots[, colorScaleName]
					} else if (!is.factor(annotVect)) {
						stop(
							colorScaleName,
							" is not factors or numeric, please check the ",
							"sample annotation table."
						)
					}
					if (sum(!levels(annotVect) %in% names(colorScale)) > 0)
						stop(
							"Levels of ",
							colorScaleName,
							" are existing in sample annotation table but ",
							"not in provided color scale."
						)
					newColorScales[[colorScaleName]] <-
						colorScale[levels(annotVect)]
				} else{
					#numeric
					if (!is.numeric(annotVect))
						stop(
							colorScaleName,
							" is not numeric but encoded as numeric ",
							"(color vector has no names)"
						)
					if (is.function(colorScale) &
							!returnContinuousFun)
						stop(
							"You must not provide function in colorScales ",
							"if returnContinuousFun=FALSE"
						)
					if (!is.function(colorScale) &
							returnContinuousFun) {
						newColorScales[[colorScaleName]] <-
							computeColorScaleFun(colorScale,
																	 values = annotVect,
																	 returnColorFun = TRUE,
																	 ...)
					} else{
						newColorScales[[colorScaleName]] <- colorScale
					}
				}
			}
			colorScalesToGen <-
				setdiff(colorScalesToGen, names(newColorScales))
		}
		cN <- 1
		cF <- 1
		for (colorScaleName in colorScalesToGen) {
			annotVect <- annots[, colorScaleName]
			if (is.numeric(annotVect)) {
				if (returnContinuousFun) {
					newColorScales[[colorScaleName]] <-
						computeColorScaleFun(
							continuousPalettes[[cN]],
							values = annotVect,
							returnColorFun = TRUE,
							...
						)
				} else{
					newColorScales[[colorScaleName]] <- continuousPalettes[[cN]]
				}
				cN <- cN + 1
				if (cN > length(continuousPalettes))
					cN <- 1
			} else{
				annots[, colorScaleName] <-
					as.factor(as.character(annots[, colorScaleName]))
				annotVect <- annots[, colorScaleName]
				newColorScales[[colorScaleName]] <-
					discreteFuns[[cF]](nlevels(annotVect))
				names(newColorScales[[colorScaleName]]) <-
					levels(annotVect)
				cF <- cF + 1
				if (cF > length(discreteFuns))
					cF <- 1
			}
		}
		newColorScales
	}

#' Add line break between factors and remove line in the middle if x axis is
#' discrete.
#'
#' @param gg ggplot object
#' @param borderSize Single numeric value. Size width of the line.
#' @param borderColor Single character value. Color of the line.
#'
#' @return A ggplot objet.
#' @export
#'
#' @examples
#' library(ggplot2)
#' g<-ggplot(data.frame(x=c("A","A","B","B","B","C")),aes(x=x))+geom_bar()
#' g
#' ggBorderedFactors(g)
#' ggBorderedFactors(g,borderColor="white",borderSize=1.5)
ggBorderedFactors<-function (gg, borderSize = 0.75, borderColor = "grey75")
{
	nX <- nlevels(as.factor(gg$data[, quo_name(gg$mapping$x)]))
	gg + geom_vline(xintercept = seq(1.5, nX - 0.5, 1), linewidth = borderSize,
									color = borderColor) + scale_x_discrete(expand = c(0,
																																		 0.5, 0, 0.5)) + theme(panel.grid.major.x = element_line(colour = NA),
																																		 											panel.grid.minor.x = element_line(colour = NA), )
}



#' Colors for a qualitative scale
#'
#' @param n Number of wanted colors
#'
#' @return A vector of colors
#' @export
#'
#' @examples
#' gsds_colors() |> scales::show_col()
#' gsds_colors(n=5) |> scales::show_col()
#' gsds_colors(n=40) |> scales::show_col()
gsds_colors <- function(n = 20) {
	myCOlors <- c(
		"#E52421",
		"#66B32E",
		"#2A4B9B",
		"#6EC6D9",
		"#F3E600",
		"#A6529A",
		"#7C1623",
		"#006633",
		"#29235C",
		"#0084BC",
		"#E6007E",
		"#F49600",
		"#E3E3E3",
		"#626F72",
		"#040505",
		"#E74B65",
		"#95B37F",
		"#683C11",
		"#F8BAA0",
		"#DD8144"
	)
	if (n <= 20) {
		return(myCOlors[seq_len(n)])
	} else{
		return(
			extendColorPalette(
				n = n,
				colors = myCOlors,
				sortColorIn = TRUE,
				sortColorOut = TRUE
			)
		)
	}
}



#' Best theoretical color palette (wrapper for qualpal)
#'
#' @inheritParams qualpalr::qualpal
#' @return Colors in hex format.
#' @export
#'
#' @examples
#' mostDistantColor(3)
#'
mostDistantColor <-
	function(n,
					 colorspace = list(h = c(0, 360), s = c(0.1, 0.9), l = c(0.1, 0.9)),
					 cvd = c(protan = 0, deutan = 0, tritan = 0)) {
		if (n == 1)
			return("#000000")
		# test if qualpalr is installed
		if (!requireNamespace("qualpalr", quietly = TRUE)) {
			warning("qualpalr is not installed, ",
							"please install for getting better colors.")
			return(ggplotColours(n))
		} else{
			qualpalr::qualpal(
				n = n,
				colorspace = colorspace,
				cvd = cvd
			)$hex
		}

	}



#' Add or remove colors to an existing palette by interpolation
#'
#' @param n Number of returned colors.
#' @param colors Vector of colors.
#' @param sortColorIn Order color vector by similarity before the interpolation.
#' @param sortColorOut Order color vector by dissimilarity after the
#'   interpolation.
#'
#' @return Vector of colors.
#' @export
#'
#' @examples
#' extendColorPalette(9,  colors=c("red","green","blue")) |> scales::show_col()
#' extendColorPalette(9,  colors=c("red","green","blue"),
#'     sortColorIn=TRUE, sortColorOut=TRUE) |> scales::show_col()
extendColorPalette <- function(n,
															 colors = c(
															 	"#E52421",
															 	"#66B32E",
															 	"#2A4B9B",
															 	"#6EC6D9",
															 	"#F3E600",
															 	"#A6529A",
															 	"#7C1623",
															 	"#006633",
															 	"#29235C",
															 	"#0084BC",
															 	"#E6007E",
															 	"#F49600",
															 	"#E3E3E3",
															 	"#626F72",
															 	"#040505",
															 	"#E74B65",
															 	"#95B37F",
															 	"#683C11",
															 	"#F8BAA0",
															 	"#DD8144"
															 ),
															 sortColorIn = FALSE,
															 sortColorOut = FALSE) {
	if (sortColorIn)
		colors <- sortColorByDistance(colors)
	colorFun <-
		circlize::colorRamp2(breaks = seq(0, 1, length.out = length(colors)),
												 colors = colors)
	colorOut <- colorFun(seq(0, 1, length.out = n))
	if (sortColorOut)
		colorOut <-
		sortColorByDistance(colorOut, byDissimilarity = TRUE)
	colorOut
}

#' Sort a vector of color by their similarity
#'
#' @param colorVector Vector of colors.
#' @param byDissimilarity Order by dissimilarity instead of similarity
#'
#' @return A vector of colors, sorted.
#' @export
#'
#' @details
#' The function proceeds as follows:
#' 1. Convert the color vector to the Lab color space.
#' 2. Compute the distance between each color.
#' 3. Order the colors using a hierarchical clustering.
#'
#' @examples
#' colors <- c("#000066","#660000","#006600","#0000FF","#FF0000","#00FF00")
#' sortColorByDistance(colors) |> scales::show_col()
#' colors <- c("#000066","#0000FF","#660000","#FF0000","#006600","#00FF00")
#' sortColorByDistance(colors, byDissimilarity=TRUE) |> scales::show_col()
sortColorByDistance <-
	function(colorVector, byDissimilarity = FALSE) {
		d <- col2rgb(colorVector) |> t()
		d <-
			grDevices::convertColor(d, from = "sRGB", to = "Lab") |>
			dist(method = "manhattan")
		if (byDissimilarity)
			d <- max(d) - d
		colorVector[hclust(d)$order]
	}




#' Create breaks for a custom ggplot scale
#'
#' @param x A vector of numeric.
#' @param scale A scale object as produced by scales::trans_new
#' @param m The number of desired breaks
#'
#' @return A vector of numeric for the major breaks
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' library(ggplot2)
#' custom_scale<-scales::trans_new(name = "invLog10",
#'     transform = function(x) -log10(x),
#'     inverse = function(x) 10^(-x), domain = c(0, Inf))
#' breaks <- ggplotBreak(DEgenesPrime_Naive$pvalue, custom_scale)
#' ggplot(DEgenesPrime_Naive, aes(x=log2FoldChange, y=pvalue)) +
#'     geom_point() +
#'     scale_y_continuous(trans = custom_scale, breaks = breaks)
ggplotBreak <- function(x, scale, m = 5) {
	transValues <- scale$transform(x)
	breaks <-
		labeling::extended(min(transValues), max(transValues), m = m)
	scale$inverse(breaks)
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
	if ((diff(h) %% 360) < 1)
		h[2] <- h[2] - 360 / n
	hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

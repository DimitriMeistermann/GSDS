whichTop<-function(x, top = 5, decreasing = TRUE) {
	order(x, decreasing = decreasing)[seq_len(top)]
}


fastWrite <- function (x,
											 fileName = "default.tsv",
											 headRow = "Name",
											 row.names = TRUE,
											 col.names = TRUE,
											 dec = ".",
											 sep = "\t",
											 ...)
{
	if (is.null(rownames(x)))
		row.names <- FALSE
	if (is.null(colnames(x)))
		col.names <- FALSE
	if (row.names) {
		x <- cbind(rownames(x), x)
		if (col.names)
			colnames(x)[1] <- headRow
	}
	data.table::fwrite(
		x = data.frame(x),
		file = fileName,
		sep = sep,
		row.names = FALSE,
		col.names = col.names,
		quote = FALSE,
		dec = dec,
		...
	)
}

checkPackage <- function(pck){
	if (!(requireNamespace(pck, quietly = TRUE))) {
		stop("Error, please download the ",pck," package. Example: BiocManager::install(\"",pck,"\")")
	}
}


#' Return ordered index of element with top n value.
#'
#' @param x Numeric vector.
#' @param top Integer. Number of element to be returned.
#' @param decreasing Logical, return top element by decreasing or increasing
#'   order.
#'
#' @return A vector of integer.
#' @export
#'
#' @examples
#' x<-c(1,5,6,10,5.2,3,8)
#' whichTop(x)
#' whichTop(x,top=3)
#' whichTop(x,decreasing=FALSE)
whichTop <- function(x, top = 5, decreasing = TRUE) {
	order(x, decreasing = decreasing)[seq_len(top)]
}



#' Draw n samples from each population and return it has a named vector
#'
#' @param sampleNames A character vector of the name of the samples.
#' @param group A factor or a character vector of the same length as
#'   `sampleNames`. Describe the population of each sample.
#' @param maxDrawSize Maximum number of observation to draw per group.
#' @param minDrawSize Minimum number of observation to draw per group.
#' @param replace Logical. Should the sample be drawn with replacement.
#'
#' @return A named vector with `sampleNames` with the population of each sample.
#' @export
#'
#' @examples
#' sampleNames<-paste0("sample",1:20)
#' group<-c(rep("A",5),rep("B",12),rep("C",3))
#' drawSamplePerGroup(sampleNames,group)
#' drawSamplePerGroup(sampleNames,group,minDrawSize=6)
#' drawSamplePerGroup(sampleNames,group,maxDrawSize=2)
#'
#' # If minDrawSize > number of sample in a group and replace=TRUE, samples
#' # will be drawn with replacement in this group
#' drawSamplePerGroup(sampleNames,group,minDrawSize=6, replace=TRUE)
drawSamplePerGroup<-function(sampleNames, group, maxDrawSize = NULL,
														 minDrawSize = NULL,replace=FALSE){
	if (is.null(group))
		return(sampleNames)
	if (is.null(names(group)))
		names(group)<-sampleNames
	if(! (is.factor(group) | is.character(group)))
		stop("group must be a factor or a character vector")
	drawSize <- min(table(group))
	if (!is.null(maxDrawSize))
		drawSize <- min(drawSize, maxDrawSize)
	if (!is.null(minDrawSize))
		drawSize <- max(drawSize, minDrawSize)
	drawCells <- unlist(lapply(unique(group), function(lvl) {
		cells <- sampleNames[group == lvl]
		if(replace & length(cells)<drawSize){
			return(sample(cells, drawSize, replace = TRUE))
		} else {
			sample(cells, min(drawSize, length(cells)))
		}
	}))
	return(group[drawCells])
}


#' Transform a range of value to another by a linear relationship.
#'
#' @description Similar to the javascript function `d3.scaleLinear()`.
#'
#' @param vals A numeric vector. Values to be transposed in the new range.
#' @param newRange A vector of two numeric values: the new minimum and maximum.
#' @param returnFunction Logical. Return the linear scale as a function instead
#'   of the transposed values in a new scale. If set to `TRUE`, `vals` argument
#'   can be also a vector of 2 numeric corresponding to the minimum and maximum
#'   of the old range.
#' @return A vector of value or a function if `returnFunction=TRUE`.
#' @export
#'
#' @examples
#' oldValues<-seq_len(10)
#' linearScale(oldValues,c(0,1),returnFunction = FALSE)
#' scaleFun<-linearScale(c(1,10),c(0,1),returnFunction = TRUE)
#' scaleFun(oldValues)
linearScale <- function(vals, newRange, returnFunction = TRUE) {
	if (!is.numeric(vals))
		stop("x should be a vector of numerics")
	if (length(newRange) != 2 |
			!is.numeric(newRange))
		stop("newRange should be a vector of 2 numerics")

	oldMin <- min(vals)
	oldMax <- max(vals)
	newMin <- newRange[1]
	newMax <- newRange[2]

	mfac <- (newMax - newMin) / (oldMax - oldMin)
	scaleFun <- function(x)
		newMin + (x - oldMin) * mfac

	if (returnFunction) {
		scaleFun
	} else{
		scaleFun(vals)
	}
}



#' Convert a list to a named factor vector
#'
#' @param listOfVector A named list. Each element must contain a character
#'   vector.
#'
#' @return A named factor vector.
#' @export
#'
#' @examples
#' VectorListToFactor(list(a=c("x1","x2"),b=c("x3","x4"),c=c("x5","x6","x7")))
#'
#' @seealso factorToVectorList
VectorListToFactor <- function(listOfVector) {
	res <-
		factor(unlist(lapply(seq_along(listOfVector), function(i)
			rep(names(listOfVector)[i], length(listOfVector[[i]])))),
			levels = names(listOfVector))
	names(res) <- unlist(listOfVector)
	res
}

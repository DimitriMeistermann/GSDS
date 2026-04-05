
hierarchicalClustering <-
	function(x,
					 transpose = TRUE,
					 method.dist = "euclidean",
					 method.hclust = "ward.D2",
					 nboot = 10,
					 PCAfirst = FALSE,
					 nDimPCA = 10) {
		if (transpose) {
			x <- t(x)
		}
		if (PCAfirst) {
			x <- fastPCA(x,
									 transpose = FALSE,
									 scale = FALSE,
									 nPC = 10
			)$x
		}
		if (is.function(method.dist)) {
			resDist <- method.dist(x)
		} else if (method.dist == "pearson") {
			resDist <- corrDist(x)
		} else {
			resDist <- dist(x, method = method.dist)
		}
		resClust <-
			stats::hclust(resDist, method = method.hclust)
		return(resClust)
	}



#' Compute correlation distance
#'
#' @param x A matrix of numeric. The function will compute the distance between
#'   the rows (same as `dist`).
#' @param method Character. Name of the method to compute correlation. See
#'   `method` argument from `cor`.
#' @param ... Other arguments to be passed to `cor`.
#'
#' @return An object of class "dist".
#' @export
#'
#' @examples
#' data("iris")
#' corrDist(t(iris[,seq_len(3)]))
corrDist <- function(x, method = "pearson", ...) {
	x <- Matrix::t(x)
	x <- stats::cor(x, method = method)
	return(as.dist( (1 - x)/2 ))
}

#' Compute cosine distance
#'
#' @param x  A matrix of numeric. The function will compute the distance between
#'   the rows (same as `dist`).
#'
#' @return An object of class "dist".
#' @export
#'
#' @examples
#' data("iris")
#' cosineDist(t(iris[,seq_len(3)]))
cosineDist<-function(x){
	sim_matrix <- coop::cosine(Matrix::t(x))
	as.dist(1 - sim_matrix)
}


#' Compute covariance distance
#'
#' @param x  A matrix of numeric. The function will compute the distance between
#'   the rows (same as `dist`).
#'
#' @return  An object of class "dist".
#' @export
#'
#' @examples
#' data("iris")
#' covDist(t(iris[,seq_len(3)]))
covDist<-function(x){
	x<-cov(
		Matrix::t(x)
	)
	max(x)-x |> as.dist()
}

#' Scale per row
#'
#' @param data A matrix or dataframe of numerics.
#' @param center either a logical value or numeric-alike vector of length equal
#'   to the number of columns of x, where ‘numeric-alike’ means that
#'   as.numeric(.) will be applied successfully if is.numeric(.) is not true.
#' @param scaled either a logical value or a numeric-alike vector of length
#'   equal to the number of columns of x.
#'
#' @return For scale.default, the centered, scaled matrix. The numeric centering
#'   and scalings used (if any) are returned as attributes "scaled:center" and
#'   "scaled:scale"
#' @export
#'
#' @examples
#' rowScale(matrix(rnorm(100),ncol = 5))
rowScale <- function(data,
										 center = TRUE,
										 scaled = FALSE) {
	data <- t(data)
	data <- t(scale(data, center = center, scale = scaled))
	return(data)
}



#' Approximation of area under ROC curve
#'
#' @description By Miron Kursa https://mbq.me
#'
#' @param score a vector of numeric representing the measure of a feature in a
#'   set of samples.
#' @param boolVect a vector of logical, same size as score. Is the sample in the
#'   target group?
#'
#' @return A numeric value. Close to 1 = perfect marker, around 0.5 = as good as
#'   random values, close to 0 = perfect anti-marker.
#' @export
#'
#' @examples
#' data(iris)
#' auroc(iris$Sepal.Length,iris$Species=="virginica")
auroc <- function(score, boolVect) {
	n1 <- sum(!boolVect)
	n2 <- sum(boolVect)
	U    <- sum(rank(score)[!boolVect]) - n1 * (n1 + 1) / 2
	return(1 - U / n1 / n2)
}

#' Quick approximation of area under ROC curve
#'
#' @description By Miron Kursa https://mbq.me
#'
#' @param score a vector of numeric representing the measure of a feature in a
#'   set of samples.
#' @param boolVect a vector of logical, same size as score. Is the sample in the
#'   target group?
#'
#' @return A numeric value. Close to 1 = perfect marker, around 0.5 = as good as
#'   random values, close to 0 = perfect anti-marker.
#' @export
#'
#' @examples
#' data(iris)
#' qauroc(iris$Sepal.Length,iris$Species=="virginica")
qauroc <- function(score, boolVect) {
	n1 <- sum(!boolVect)
	n2 <- sum(boolVect)
	aurocCPP(as.numeric(score),
					 as.logical(boolVect),
					 as.integer(n1),
					 as.integer(n2))
}

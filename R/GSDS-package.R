#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import grid
#' @import ggplot2
#' @import ComplexHeatmap
#' @importFrom SummarizedExperiment
#'  SummarizedExperiment
#'  assay
#'  colData
#'  rowData
#'  rowData<-
#'  assays
#' @importFrom AnnotationDbi keys
#' @importFrom AnnotationDbi select
#' @importFrom fgsea fgseaMultilevel
#' @importFrom fgsea reactomePathways
#' @importFrom stats
#'  cor
#'  cov
#'  dist
#'  as.dist
#'  lm
#'  p.adjust
#'  sd
#'  na.omit
#'  as.formula
#'  setNames
#'  hclust
#'  .lm.fit
#'  pbinom
#'  phyper
#'  quantile
#' @importFrom grDevices col2rgb hcl
#' @importFrom RSQLite dbConnect SQLite dbGetQuery
#' @importFrom utils data
#' @importFrom S4Vectors DataFrame
#' @useDynLib GSDS, .registration=TRUE
## usethis namespace: end
NULL


#' species_df
#' @name species_df
#' @docType data
#' @author Dimitri Meistermann
#' @keywords data
#' @description
#' A data frame describing available species (available Organism-level Databases),
#' Columns: package, species, simple_species, kegg_code, id_type, tax_id"
#' @format A data frame containing 21 rows for 6 columns
NULL

#' keggDB
#' @name keggHuman
#' @docType data
#' @author Kanehisa M
#' @keywords data
#' @references \url{https://doi.org/10.1093/nar/28.1.27}
#' @description
#' KEGG gene sets for human
#' @format A list where first element is a list of 369 gene sets
NULL


#' sampleAnnot
#'
#' @name sampleAnnot
#' @docType data
#' @author Dimitri Meistermann
#' @references \url{https://doi.org/10.1038/s41467-017-02107-w}
#' @keywords data
#' @description Sample Annotation sheet (often referred as colData) for 35 bulk
#' RNA-Seq samples used in Kilens, Meistermann et al. 2018.
#'
#'
#' @format A data.frame containing 35 observations for 18 features.
NULL

#' bulkLogCounts
#'
#' @name bulkLogCounts
#' @docType data
#' @author Dimitri Meistermann
#' @references \url{https://doi.org/10.1038/s41467-017-02107-w}
#' @keywords data
#' @description Log-expression count table for 35 samples sequenced from a bulk
#' 3'digital RNA-Seq protocol used in Kilens, Meistermann et al. 2018.
#' @format A matrix containing 35 observations for 16959 features/genes
NULL


#' DEgenesPrime_Naive
#'
#' @name DEgenesPrime_Naive
#' @docType data
#' @author Dimitri Meistermann
#' @references \url{https://doi.org/10.1038/s41467-017-02107-w}
#' @keywords data
#' @description
#' Results of the differential expression analysis between primed and naive
#' human pluripotent stem cells in Kilens, Meistermann et al. 2018.
#' @format A data.frame containing 16959 observations for 7 features
NULL

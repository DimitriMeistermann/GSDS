
#' Return gene id correspondence, GO species code and KEGG species code
#'
#' @param sample.species Character. Shortname of the species as described in `data("bods")`.
#'
#' @return A list describing specific data for the species (gene IDs, annotation package...).
#' @export
#'
#' @examples
#' require(BiocManager)
#' require(org.Hs.eg.db)
#'
#' getSpeciesData("Human") |> str()
#'
#' #Available species
#' data("bods", package = "gage")
#' bods[,"species"]
getSpeciesData <-
    function(sample.species = "Human") {
        bods <- ""
        data("bods", package = "gage", envir = environment())
        species <- list()
        species.data <- data.frame(bods)
        species.index <- which(species.data$species == sample.species)
        if (length(species.index) != 1)
            stop(
                "Wrong species name, type \ndata(bods, package = 'gage')\nbods[,\"species\"]\nto see available species"
            )
        species$package <-
            as.character(species.data[species.index, "package"])
        species$kegg <-
            as.character(species.data[species.index, "kegg.code"])
        species$go <-
            strsplit(as.character(species$package),
                     split = ".",
                     fixed = TRUE)[[1]][2]
        if (!(requireNamespace(species$package, quietly = TRUE))) {
            stop(paste0("Error, please download the following annotation package, then retry: ",
                        species$package ))
        }
        getOfSpeciesPackage<-AnnotationDbi::get(species$package,envir=getNamespace(species$package))
        species$GeneIdTable <-
            AnnotationDbi::select(
                getOfSpeciesPackage,
                keys = AnnotationDbi::keys(getOfSpeciesPackage, "ENTREZID") ,
                columns = c("ENTREZID", "SYMBOL", "ENSEMBL")
            ) |> suppressMessages()
        species$species <- sample.species
        return(species)
    }



#' Download database of term/pathway for enrichment analyses.
#'
#' @param geneSym A vector of gene symbols.
#' @param geneEntrez NULL or a vector of gene Entrez ID.
#' @param idGeneDF Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom.
#' @param customAnnot  Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param returnGenesSymbol Logical, return gene symbol in each term instead of Entrez ID.
#' @param filter_genes Logical, filter genes that are not in the database.
#'
#' @return
#' A list where each element is a database that contain a list of term with the associated gene symbols.
#' @export
#' @import AnnotationDbi
#' @examples
#' data("bulkLogCounts", package="oob")
#' enrichDBs<-getDBterms(rownames(bulkLogCounts),species="Human",database=c("kegg","goBP"))
getDBterms<-function(geneSym = NULL,geneEntrez=NULL, idGeneDF=NULL, speciesData=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),customAnnot=NULL,
                     keggDisease=FALSE,species="Human",returnGenesSymbol=TRUE, filter_genes = TRUE){
    select<-AnnotationDbi::select
    validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
    if(sum(database%in%validDBs)<length(database)) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
    if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")
    if(is.null(speciesData)){
        speciesData<-getSpeciesData(species)
    }else{
        species<-speciesData$species
    }
    if(is.null(idGeneDF)) idGeneDF<-speciesData$GeneIdTable
    options(warn=-1)
    if(is.null(geneEntrez)){
        if(is.null(geneSym)) stop("You must give a value to geneSym or geneEntrez")
        geneEntrez<-oob::ConvertKey(geneSym,tabKey = idGeneDF,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
        geneEntrez<-geneEntrez[!is.na(geneEntrez)]
    }
    db_terms<-list()
    if(is.list(customAnnot)){
        db_terms$custom<-lapply(customAnnot,function(x){
            new_x<-oob::ConvertKey(x,tabKey = idGeneDF,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
            new_x[!is.na(new_x)]
        })
    }
    if(!(length(database)<=1 & database[1]=="custom")){
        if("reactom"%in%database){
            if (!(requireNamespace("reactome.db", quietly = TRUE))) {
                stop("Error, please download 'reactome.db' package, then retry")
            }
            db_terms$reactom<- fgsea::reactomePathways(geneEntrez)
            db_terms$reactom<-db_terms$reactom[unique(names(db_terms$reactom))]
        }
        if("kegg"%in%database){
            kg.species <- gage::kegg.gsets(speciesData$kegg, id.type="entrez")
            db_terms$kegg<- if(keggDisease) kg.species$kg.sets else kg.species$kg.sets[kg.species$sigmet.idx]
        }
        if("go"%in%substr(database,1,2)){
            go.species <- gage::go.gsets(tolower(species), id.type="entrez")
            if("goBP"%in%database) db_terms$goBP<-go.species$go.sets[go.species$go.subs$BP]
            if("goMF"%in%database) db_terms$goMF<-go.species$go.sets[go.species$go.subs$MF]
            if("goCC"%in%database) db_terms$goCC<-go.species$go.sets[go.species$go.subs$CC]
        }
    }
    options(warn=0)

    if(filter_genes){
        db_terms<-lapply(db_terms,function(db) lapply(db,function(x) x[x%in%geneEntrez]))
    }

    if(returnGenesSymbol){
        lapply(db_terms,function(db) lapply(db,oob::ConvertKey,tabKey=idGeneDF,colOldKey = "ENTREZID",colNewKey = "SYMBOL"))
    }else{
        db_terms
    }
}




#' Export enrichment results with the gene list associated to each term/row.
#'
#' @param enrichResults Dataframe, usually from Enrich functions (for example `enrich.ora`) with `returnGenes=TRUE`.
#' @param file Path to file to write.
#' @param sep Field separator.
#' @param col.names Either a logical value indicating whether the column names of x are to be written along with x, or a character vector of column names to be written. See the section on ‘CSV files’ for the meaning of col.names = NA.
#' @param row.names Either a logical value indicating whether the row names of x are to be written along with x, or a character vector of row names to be written.
#' @param geneCol The column that contain the list of gene of the term.
#' @param ... Other parameters passed to write.table.
#'
#' @return Write a text file.
#' @export
#'
exportEnrich <-
    function(enrichResults,
             file,
             sep = "\t",
             col.names = TRUE,
             row.names = FALSE,
             geneCol = "genes",
             ...) {
        genesPerRow <- enrichResults[, geneCol]
        if (!is.list(genesPerRow))
            stop(geneCol, " must be a list")
        if (!all(vapply(genesPerRow, is.character, FUN.VALUE = logical(1))))
            stop(geneCol, " must be a list of character vectors")


        enrichResults[, geneCol] <- NULL

        enrichResults[,geneCol] <-
            vapply(genesPerRow, FUN.VALUE = character(1), function(x) {
                return(paste0(x, collapse = sep))
            })

        oob::fastWrite(
            enrichResults,
            file,
            sep = sep,
            col.names = col.names,
            row.names = row.names,
            ...
        )
    }


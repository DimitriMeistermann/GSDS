
#' Return gene id correspondence, GO species code and KEGG species code
#'
#' @param sample.species Character. Shortname of the species as described in `data("bods")`.
#' @param updateSpeciesPackage Logical. Download or update automatically the annotation org.db package corresponding to the species.
#'
#' @return A list describing specific data for the species (gene IDs, annotation package...).
#' @export
#'
#' @examples
#' library(gage)
#' data(bods)
#' bods
#' getSpeciesData("Human")
#' getSpeciesData("Mouse")
getSpeciesData<-function(sample.species="Human",updateSpeciesPackage=FALSE){
    data("bods",package = "gage")
    species<-list()
    species.data<-data.frame(bods)
    species.index<-which(species.data$species==sample.species)
    if(length(species.index)!=1) stop("Wrong species name, type \ndata(bods)\nbods[,\"species\"]\nto see available species")
    species$package<-as.character(species.data[species.index,"package"])
    species$kegg<-as.character(species.data[species.index,"kegg.code"])
    species$go<-strsplit(as.character(species$package),split = ".",fixed = TRUE)[[1]][2]
    if(updateSpeciesPackage | !(require(species$package,character.only = TRUE))){
        print(paste0("Downloading species package: ",species.data$package))
        BiocManager::install(species$package, update=FALSE)
    }
    require(species$package,character.only = TRUE)
    species$GeneIdTable<-AnnotationDbi::select(get(species$package),keys = AnnotationDbi::keys(get(species$package),"ENTREZID") , columns = c("ENTREZID","SYMBOL","ENSEMBL")) |> suppressMessages()
    species$species<-sample.species
    return(species)
}




#' Download database of term/pathway for enrichment analyses.
#'
#' @param geneSym A vector of gene symbols.
#' @param geneEntrez NULL or a vector of gene Entrez ID.
#' @param corrIdGenes Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom.
#' @param customAnnot  Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param returnGenesSymbol Logical, return gene symbol in each term instead of Entrez ID.
#'
#' @return
#' A list where each element is a database that contain a list of term with the associated gene symbols.
#' @export
#' @import AnnotationDbi
#' @examples
#' data("bulkLogCounts")
#' enrichDBs<-getDBterms(rownames(bulkLogCounts),species="Human",database=c("kegg","reactom"))
getDBterms<-function(geneSym,geneEntrez=NULL, corrIdGenes=NULL, speciesData=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),customAnnot=NULL,
                     keggDisease=FALSE,species="Human",returnGenesSymbol=TRUE){
    select<-AnnotationDbi::select
    validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
    if(sum(database%in%validDBs)<length(database)) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
    if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")
    if(is.null(speciesData)){
        speciesData<-getSpeciesData(species)
    }else{
        species<-speciesData$species
    }
    if(is.null(corrIdGenes)) corrIdGenes<-speciesData$GeneIdTable
    options(warn=-1)
    if(is.null(geneEntrez)){
        geneEntrez<-ConvertKey(geneSym,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
        geneEntrez<-geneEntrez[!is.na(geneEntrez)]
    }
    db_terms<-list()
    if(is.list(customAnnot)){
        db_terms$custom<-lapply(customAnnot,function(x){
            new_x<-ConvertKey(x,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
            new_x[!is.na(new_x)]
        })
    }
    if(!(length(database)<=1 & database[1]=="custom")){
        if("reactom"%in%database){
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

    if(returnGenesSymbol){
        lapply(db_terms,function(db) lapply(db,ConvertKey,tabKey=corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL"))
    }else{
        db_terms
    }
}


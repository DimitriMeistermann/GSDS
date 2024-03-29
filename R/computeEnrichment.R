#' Functional class scoring enrichment (fgsea algorithm)
#'
#' @param x vector or dataframe/matrix of one column. Values are used to classify the genes, example, it can be Log2(Fold-Change0). Genes are contained in the names/rownames of the vector/dataframe. Example of valid x: x<-rnorm(n = 4); names(x)<-c("GATA2","SOX17","KLF4","POU5F1")
#' @param idGeneDF Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom
#' @param maxSize maximum number of gene in each term
#' @param minSize Minimum number of gene in each term
#' @param customAnnot custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param returnGenes  return genes that were the most important for the enrichment of term
#' @param keggDisease Logical. Retain kegg disease term ?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param db_terms A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
#' If this argument is not NULL, no additional database are downloaded.
#' @param speciesData  object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#' @param ... Additional parameters that are passed to fgsea
#'
#' @return
#' A data frame with the following columns:
#' - pathway: name of the pathway/term
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#' - log2err: the expected error for the standard deviation of the P-value logarithm
#' - ES: enrichment score, same as in Broad GSEA implementation
#' - NES: enrichment score normalized to mean enrichment of random samples of the same size
#' - size: number of gene in the term after removing genes not present
#' - genes (if `returnGenes`). Vector of genes of the term.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive",package = "oob")
#' fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),
#'     DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
#' resEnrich<-enrich.fcs(fcsScore,database = "kegg",species = "Human")
#' head(resEnrich)
enrich.fcs<-function(x, idGeneDF=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
                     maxSize=500,minSize=2,customAnnot=NULL,returnGenes=FALSE,
                     keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL,...){
    if(is.data.frame(x) | is.matrix(x)){
        if(is.null(rownames(x))) stop("If input format is a matrix/dataframe, row should be named with gene symbols")
        tempx<-x
        x<-tempx[,1]
        names(x)<-rownames(tempx)
    }

    if(is.null(names(x))) stop("Values should be named with gene symbols")
    if(!is.numeric(x)) stop("Values must be numeric")
    if(is.null(db_terms)) db_terms<-getDBterms(geneSym=names(x), idGeneDF=idGeneDF,database=database,
                                               customAnnot=customAnnot,keggDisease=keggDisease,species=species)
    if(length(db_terms)==0) stop("Error, no term in any database was found")
    res<-list()
    for(db in names(db_terms)){
        res[[db]]<-suppressWarnings(fgsea::fgseaMultilevel (db_terms[[db]], x ,minSize=minSize,maxSize=maxSize,eps = 0,...))
        res[[db]]<-res[[db]][order(res[[db]]$padj),]
        res[[db]]$database<-db
        res[[db]]$leadingEdge<-NULL
        if(returnGenes) res[[db]]$genes <- db_terms[[db]][res[[db]]$pathway]
    }
    res<-do.call("rbind", res)
    res$padj<-p.adjust(res$pval,method = "BH")
    return(data.frame(res))
}


#' Over Representation Analysis (enrichment, Fischer tests)
#'
#' @param x vector or dataframe/matrix of one column.
#' Values are booleans and say if gene is from the list of interest or not.
#' Genes are contained in the names/rownames of the vector/dataframe.
#' Example of valid x: x<-c(TRUE,TRUE,FALSE,FALSE); names(x)<-c("GATA2","SOX17","KLF4","POU5F1").
#' In this case, GATA2, SOX17, KLF4, POU5F1 are the universe of gene and GATA2 and SOX17 are the genes of interest
#' @param idGeneDF  Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom
#' @param minSize Minimum number of gene in each term.
#' @param maxSize Maximum number of gene in each term.
#' @param returnGenes Return genes that were the most important for the enrichment of term.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param customAnnot Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param db_terms A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#'
#' @return
#' A dataframe with the following columns:
#' - term: name of the term
#' - pval: an enrichment p-value
#' - OD: Odds Ratio of the enrichment
#' - padj: a BH-adjusted p-value
#' - nGeneOfInterest: number of gene of interest in the term.
#' - nGene: number of gene in the term after removing genes not present.
#' - genes (if `returnGenes`). Vector of genes of the term.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive" ,package = "oob")
#' vectorIsDE<-DEgenesPrime_Naive$isDE!="NONE";names(vectorIsDE)<-rownames(DEgenesPrime_Naive)
#' resEnrich<-enrich.ora(vectorIsDE,database = "kegg",species = "Human")
#' head(resEnrich)
enrich.ora<-function(x, idGeneDF=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
                     minSize=2,maxSize=500,returnGenes=FALSE, keggDisease=FALSE,species="Human",
                     customAnnot=NULL,db_terms=NULL,speciesData=NULL){
    validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
    if(sum(database%in%validDBs)==0) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
    if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")

    if(is.data.frame(x) | is.matrix(x)){
        tempx<-x
        x<-tempx[,1]
        names(x)<-rownames(tempx)
    }

    if(!is.logical(x)) stop("Values must be logical (TRUE or FALSE)")
    if(!is.character(names(x))) stop("Values must be named with genes symbol")

    if(is.null(db_terms)) db_terms<-getDBterms(geneSym=names(x), idGeneDF=idGeneDF,
                                               database=database,customAnnot=customAnnot,keggDisease=keggDisease,species=species)

    nInterest<-length(which(x))
    nuniverse<-length(x)

    results<-list()
    for(db in names(db_terms)){
        len_term<-sapply(db_terms[[db]],length)
        db_terms[[db]]<-db_terms[[db]][len_term>=minSize & len_term<=maxSize]
        terms<-db_terms[[db]]

        nGeneByterm<-sapply(terms,length)
        nGeneOfInterestByterm<-sapply( terms,function(term){
            return(length(which(x[term])))
        })


        results[[db]]<-data.frame(row.names = names(terms))
        results[[db]]$term <- names(terms)
        results[[db]]$database<-db

        parameterList4Enrich<-vector(mode="list", length = length(terms))
        for(i in seq_along(terms)){
            parameterList4Enrich[[i]]<-list(intersectionSize=nGeneOfInterestByterm[i], setSizes=c(nInterest,nGeneByterm[i]), universeSize=nuniverse)
        }

        resEnrich<-lapply(parameterList4Enrich,function(params) do.call("enrichSetIntersection",params))
        results[[db]]$nGene<-nGeneByterm

        results[[db]]$obsOverlap<-nGeneOfInterestByterm
        results[[db]]$expectOverlap<- sapply(resEnrich, function(x) x$expected)
        results[[db]]$OEdeviation<- sapply(resEnrich, function(x) x$OEdeviation)
        results[[db]]$pval<- sapply(resEnrich, function(x) x$pval)
        results[[db]]$padj<-0

        if(returnGenes){
            results[[db]]$genes<- db_terms[[db]]
        }
    }
    results<-do.call("rbind", results)
    results$padj<-p.adjust(results$pval,method = "BH")
    return(results)
}



#' Gene Set Differential Scoring (GSDS)
#'
#' @param geneSetActivScore A list of database with an element "eigen" containing the matrix of gene set activation score (see what `computeActivationScore` returns).
#' If NULL, this is computed automatically from the `exprMatrix` and the gene set database given via `db_terms` or requested via `database`.
#' @param exprMatrix An expression matrix (normalized log2(x+1) counts). Genes as rows and sample as columns. If `db_terms` is not given, must be named by gene symbols.
#' @param colData An annotation dataframe. Each column is a feature, each row a sample. Same number of samples than in `exprMatrix`.
#' @param contrast A vector of 3 character.
#' 1. Name of the experimental variable that have to be used for differential activation. Must be a column name of `colData`.
#' 2. Condition considered as the reference.
#' 3. Condition considered as the target group.
#' @param idGeneDF Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom.
#' @param maxSize Maximum number of gene in each term.
#' @param minSize Minimum number of gene in each term.
#' @param customAnnot  Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param db_terms A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#'
#' @return
#' A dataframe with the following columns:
#' - term: name of the term/gene set
#' - baseMean: mean of activation score in the gene set
#' - sd: standard deviation  in the gene set
#' - log2FoldChange: Log(Log Fold Change) of activation score between the two tested groups.
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#' - database: origin of the gene set
#' - size: number of gene in the term after removing genes not present.
#' @export
#'
#' @examples
#' data("bulkLogCounts", package = "oob")
#' data("sampleAnnot", package = "oob")
#'
#' keggDB<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,db_terms = keggDB)
#' resGSDS<-GSDS(geneSetActivScore = geneSetActivScore,colData = sampleAnnot,
#'     contrast = c("culture_media","T2iLGO","KSR+FGF2"),db_terms =  keggDB)
#'
#' # or
#' resGSDS<-GSDS(exprMatrix = bulkLogCounts,colData = sampleAnnot,
#'     contrast = c("culture_media","T2iLGO","KSR+FGF2"),database = "kegg")
#'
#' head(resGSDS)
GSDS <-
    function(geneSetActivScore = NULL,
             exprMatrix = NULL,
             colData,
             contrast ,
             idGeneDF = NULL,
             database = c("kegg", "reactom", "goBP", "goCC", "goMF"),
             maxSize = 500,
             minSize = 2,
             customAnnot = NULL,
             keggDisease = FALSE,
             species = "Human",
             db_terms = NULL,
             speciesData = NULL) {

        if (is.null(geneSetActivScore) &
            is.null(exprMatrix))
            stop("At least exprMatrix or geneSetEigens must be given")


        if (is.null(geneSetActivScore)) {
            if (is.null(db_terms)) {
                db_terms <-
                    getDBterms(
                        geneSym = rownames(exprMatrix),
                        idGeneDF = idGeneDF,
                        database = database,
                        customAnnot = customAnnot,
                        keggDisease = keggDisease,
                        species = species,
                        returnGenesSymbol = TRUE
                    )
            }
            geneSetActivScore <-
                computeActivationScore(exprMatrix = exprMatrix, db_terms = db_terms)

        }

        if (!is.list(geneSetActivScore) |
            (is.list(geneSetActivScore) &
             setequal(
                 c("activScoreMat", "contributionList"),
                 names(geneSetActivScore)
             ))) {
            geneSetActivScore <- list("db" = geneSetActivScore)
        }

        res <- list()

        for (db in names(db_terms)) {
            if (is.list(geneSetActivScore[[db]])) {
                activScorePerPathway <- geneSetActivScore[[db]]$activScoreMat
            } else{
                activScorePerPathway <- geneSetActivScore[[db]]
            }

            if (!setequal(rownames(colData),
                          rownames(geneSetActivScore[[1]]$activScoreMat))) {
                stop("colData and geneSetActivScore must have the same samples names, in same order")
            }

            db_terms[[db]] <-
                db_terms[[db]][colnames(activScorePerPathway)]

            res[[db]] <-
                dfres <-
                data.frame(
                    term = names(db_terms[[db]]),
                    multiLinearModel(t(activScorePerPathway), colData, contrast = contrast),
                    database = db,
                    size = sapply(db_terms[[db]], length),
                    sd = apply(activScorePerPathway, 2, sd),
                    row.names = NULL
                )
        }
        do.call("rbind", res)

    }



#' Test a linear model on each gene following an experimental design.
#'
#' @param exprData A matrix of numeric with rows as features (in the RNA-Seq context, log counts).
#' @param colData A dataframe of feature (numeric or factor) with rows as samples. Must have the same number of samples than exprData
#' @param contrast A vector of 3 character.
#' 1. Name of the experimental variable that have to be used for differential activation. Must be a column name of `colData`.
#' 2. Condition considered as the reference.
#' 3. Condition considered as the target group.
#'
#' @return
#' A dataframe with the following columns:
#' - baseMean: mean
#' - log2FoldChange: Log(Log Fold Change)  between the two tested groups.
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts", package = "oob")
#' data("sampleAnnot", package = "oob")
#' res<-multiLinearModel(bulkLogCounts,colData = sampleAnnot,
#'     contrast = c("culture_media","T2iLGO","KSR+FGF2"))
multiLinearModel <-
    function(exprData,
             colData ,
             contrast) {
        samples <- colData[, contrast[1]] %in% contrast[2:3]
        data <- exprData[, samples, drop = FALSE]
        groups <-
            droplevels(as.factor(colData[samples, contrast[1]]))
        logicGroup <- rep(FALSE, length(groups))
        logicGroup[groups == contrast[2]] <- TRUE
        regTabList <- apply(data, 1, function(x) {
            data.frame(data = x, group = logicGroup)
        })
        resList <- lapply(regTabList, function(regTab) {
            summary(lm(data ~ group, data = regTab))$coefficients[2, c(1, 4)]
        })
        res <-
            data.frame(do.call("rbind", resList))
        colnames(res) <- c("log2FoldChange", "pval")
        res <-
            cbind(data.frame(baseMean = apply(exprData[, samples], 1, mean)), res)
        res$padj <- p.adjust(res$pval, method = "BH")
        return(res)
    }

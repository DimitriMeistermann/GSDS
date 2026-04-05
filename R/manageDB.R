#' Return gene id correspondence, GO species code and KEGG species code
#'
#' @param species Character. Shortname of the species (rownames in `data("species_df")`).
#' @param species_df Eventually an updated `species_df` dataframe.
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
#' data("species_df")
#' species_df[,"species"]
getSpeciesData <-
	function(species = "Human",
					 species_df = NULL) {
		if(is.null(species_df)){
			data("species_df", envir = environment())
		} else {
			spDfGoodcol <- c('package','species','kegg_code','id_type','tax_id')
			if(!all(spDfGoodcol %in% colnames(species_df))){
				stop("species_df must contain the following columns: ",paste0(spDfGoodcol, collapse = ", "))
			}
		}
		species <- checkSpecies(species, species_df = species_df)
		speciesDat<-as.list(species_df[species,])
		speciesDat$species <- species

		checkPackage(speciesDat$package) # package is installed

		# method: get species package object
		speciesDat$getSpeciesPackage <- function(){
			AnnotationDbi::get(speciesDat$package, envir =
				getNamespace(speciesDat$package))
		}
		# method: convert gene ids
		speciesDat$convertID <- function(genes,oldID,newID,multiVals="FIRST"){
			checkIDs(c(oldID,newID),speciesDat$package)
			AnnotationDbi::mapIds(speciesDat$getSpeciesPackage(),
						 keys = genes,
						 column = newID,
						 keytype = oldID,
						 multiVals = "first") # Takes the first match if multiple exist
		}

		return(speciesDat)
	}

#' Download database of term/pathway for enrichment analyses.
#'
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom.
#' @param id_type Character. Gene ID type that has to be returned, typically values can be "ENSEMBL", "ENTREZID" or "SYMBOL", see columns(speciesData$getSpeciesPackage()) to see available IDs
#' @param GOALL Logical. If `TRUE` the GO term correspond to propagated annotations for all parents term, which make the database huge an less specific, but more complete.
#' @param species Character. Shortname of the species as described in `data("species_df")`.
#'
#' @return
#' A list where each element is a database that contain a list of term with the associated gene symbols.
#' @export
#' @import AnnotationDbi
#' @examples
#' data("bulkLogCounts")
#' enrichDBs<-getDBsets(species="Human",database=c("kegg","goBP"))
getDBsets <- function(species = "Human",
											database = c("kegg", "reactom", "goBP", "goCC", "goMF"),
											id_type = "SYMBOL",  GOALL = FALSE,
											speciesData = NULL) {
	validDBs <- c("kegg", "reactom", "goBP", "goCC", "goMF")
	if (sum(database %in% validDBs) < length(database))
		stop(paste0(
			"Error, valid values for database are: ",
			paste0(validDBs, collapse = ", ")
		))
	if (is.null(speciesData)) {
		speciesData <- getSpeciesData(species)
	}
	message(
		"Downloading ",species," gene sets from ",
		paste0(database, collapse = ", "), " returning ", id_type, " IDs"
	)
	options(warn = -1)
	DBsets <- list()

	if ("reactom" %in% database) {
		DBsets$reactom <- reactomePathways(speciesData, id_type = id_type)
	}
	if ("kegg" %in% database) {
		DBsets$kegg <- keggPathways(speciesData, id_type =)
	}
	go2dl <- grep("^go", database, value = TRUE)
	if (length(go2dl) > 0) {
		DBsets <- c(DBsets,
									goSets(speciesData, goDB = go2dl, id_type = id_type, GOALL = GOALL))
	}
	options(warn = 0)

	return(DBsets)
}


#' Export enrichment results with the gene list associated to each term/row.
#'
#' @param enrichResults Dataframe, usually from Enrich functions (for example
#'   `enrich.ora`) with `returnGenes=TRUE`.
#' @param file Path to file to write.
#' @param sep Field separator.
#' @param col.names Either a logical value indicating whether the column names
#'   of x are to be written along with x, or a character vector of column names
#'   to be written. See the section on ‘CSV files’ for the meaning of col.names
#'   = NA.
#' @param row.names Either a logical value indicating whether the row names of x
#'   are to be written along with x, or a character vector of row names to be
#'   written.
#' @param geneCol The column that contain the list of gene of the term.
#' @param ... Other parameters passed to write.table.
#'
#' @return Write a text file.
#' @export
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

        fastWrite(
            enrichResults,
            file,
            sep = sep,
            col.names = col.names,
            row.names = row.names,
            ...
        )
    }

#' Download GO term for a particular Ensembl dataset, useful for non model species
#'
#' @param ensemblDataset name (species) of the Ensembl dataset
#' @param returnGeneSymbol Return gene symbol if `TRUE` or Ensembl IDs
#'
#' @returns A list where each element is a database that contain a list of term with the associated gene symbols.
#' @export
#'
#' @examples
#' \dontrun{
#' dbGO_Mzebra <- downloadGODB_BioMart("mzebra_gene_ensembl")
#' }
downloadGODB_BioMart <- function(ensemblDataset="mzebra_gene_ensembl", returnGeneSymbol = FALSE) {
	checkPackage("biomaRt")
	ensembl <- biomaRt::useMart("ensembl", dataset = ensemblDataset) # Replace "mzebra_gene_ensembl" if the dataset name is different.
	idcol <- ifelse(returnGeneSymbol, "external_gene_name","ensembl_gene_id")
	go_data <- biomaRt::getBM(
		attributes = c(idcol, "go_id", "name_1006", "namespace_1003"), #, "description", ""ensembl_gene_id"),
		mart = ensembl
	)
	go_data <- go_data[go_data[,idcol] != "" & go_data$go_id != "",]
	go_data$go_term <- paste(go_data$go_id,go_data$name_1006,sep = " ")

	db <- list(
		"goCC" = go_data[go_data$namespace_1003 == "cellular_component",c(idcol,"go_term")],
		"goMF" = go_data[go_data$namespace_1003 == "molecular_function",c(idcol,"go_term")],
		"goBP" = go_data[go_data$namespace_1003 == "biological_process",c(idcol,"go_term")]
	)
	lapply(db, function(x) {
		return(split(x[,idcol],x$go_term))
	})
}


checkSpecies <- function(species, species_df) {
	matches <- grep(species, rownames(species_df), ignore.case = TRUE, value = TRUE)

	if (length(matches) == 0) {
		stop(sprintf("No species found matching '%s'. Check species_df for valid names.", species))
	}

	if (length(matches) > 1) {
		stop(sprintf("Multiple species match '%s': (%s). Please be more specific.",
								 species, paste(matches, collapse = ", ")))
	}

	# Return the unique full name for downstream use
	return(matches)
}

checkIsSpeciesData <- function(speciesData) {
	required <- c("package", "kegg", "go", "GeneIdTable")
	if (!all(required %in% names(speciesData))) {
		stop("speciesData is missing required components. Re-run getSpeciesData().")
	}
}

keggPathways <- function(speciesData, id_type = "SYMBOL") {
	checkPackage("KEGGREST")
	species<-speciesData$kegg_code
	pathway_list <- KEGGREST::keggList("pathway", species)
	pathway_gene_link <- KEGGREST::keggLink(species, "pathway")

	# Remove "path:" prefix from pathway IDs
	path_ids <- gsub("path:", "", names(pathway_gene_link))
	# Remove "species:" prefix (e.g., "hsa:") from gene IDs
	gene_ids <- gsub(paste0(species, ":"), "", pathway_gene_link)

	if(id_type!="ENTREZID"){
		gene_ids<-speciesData$convertID(gene_ids,"ENTREZID",id_type)
	}

	# Clean up descriptions: remove " - Homo sapiens (human)" or similar suffixes
	clean_descriptions <- gsub(paste0(" - .*"), "", pathway_list)

	# Create a lookup vector: Name = clean ID (hsa...), Value = Description
	names(clean_descriptions) <- gsub("path:", "", names(clean_descriptions))

	# Map descriptions to the current list of pathway IDs
	current_descriptions <- clean_descriptions[path_ids]

	full_pathway_names <- paste(path_ids, current_descriptions)
	gene_set_list <- split(as.vector(na.omit(gene_ids)), full_pathway_names)

	return(gene_set_list)
}

reactomePathways <- function(speciesData, id_type = "SYMBOL"){
	checkPackage("reactome.db")
	pathways <- AnnotationDbi::select(
			reactome.db::reactome.db,
			keys = AnnotationDbi::keys(speciesData$getSpeciesPackage(), "ENTREZID"),
			c("PATHID"),
			keytype = "ENTREZID"
	)
	colnames(pathways)[1] <- "GENE"

	if(id_type!="ENTREZID"){
		pathways$GENE<-speciesData$convertID(pathways$GENE,"ENTREZID",id_type)
	}

	pathways <- na.omit(pathways)

	pathways <- split(pathways$GENE, pathways$PATHID)
	pathway2name <- as.data.frame(AnnotationDbi::select(
		reactome.db::reactome.db,
		names(pathways),
		c("PATHNAME"),
		"PATHID"
	))
	pathway2name <- pathway2name[!duplicated(pathway2name$PATHID),]
	pathway2name$PATHNAME <- sub("^[^:]*: ", "", pathway2name$PATHNAME)
	pathway2name <- setNames(pathway2name[["PATHNAME"]], pathway2name[["PATHID"]])
	names(pathways) <- pathway2name[names(pathways)]

	pathways
}

goSets <- function(speciesData, goDB = c("goBP","goCC","goMF"), id_type = "SYMBOL", GOALL = FALSE){

	pkg_obj <- speciesData$getSpeciesPackage()
	db_conn <- AnnotationDbi::dbconn(pkg_obj)

	# 1. Map id_type to the correct internal SQLite table and column
	# SYMBOL is in 'gene_info', ENSEMBL is in 'ensembl'
	map_info <- switch(id_type,
										 "SYMBOL"   = list(table = "gene_info", col = "symbol"),
										 "ENSEMBL"  = list(table = "ensembl",   col = "ensembl_id"),
										 "ENTREZID" = list(table = "genes",     col = "gene_id"),
										 stop("Unsupported id_type. Use 'SYMBOL', 'ENSEMBL', or 'ENTREZID'.")
	)

	# 2. Determine GO table
	go_table <- if(GOALL) "go_all" else "go"

	# 3. SQL Join using the correct table mapping
	query <- paste0(
		"SELECT DISTINCT g.go_id, g.ontology, i.", map_info$col, " AS GENE ",
		"FROM ", go_table, " AS g ",
		"JOIN ", map_info$table, " AS i ON g._id = i._id"
	)

	all_data <- RSQLite::dbGetQuery(db_conn, query)

	# 4. Get Term Names from GO.db
	go_db_conn <- AnnotationDbi::dbconn(GO.db::GO.db)
	term_info <- RSQLite::dbGetQuery(go_db_conn, "SELECT go_id, term FROM go_term")
	term_map <- setNames(paste(term_info$go_id, term_info$term), term_info$go_id)

	# 5. Format the output
	res <- lapply(goDB, function(db){
		current_onto <- substr(db, 3, 4)
		subset_df <- all_data[all_data$ontology == current_onto, ]

		# Ensure uniqueness in R as a final safety check
		subset_df <- unique(subset_df[, c("go_id", "GENE")])

		split_list <- split(subset_df$GENE, subset_df$go_id)
		names(split_list) <- term_map[names(split_list)]

		return(split_list[!is.na(names(split_list))])
	})

	names(res) <- goDB
	return(res)
}

checkIDs<-function(IDs,package){
	getpckg <- AnnotationDbi::get(package, envir =
										 	getNamespace(package))
	lapply(IDs, function(x){
		if(!x %in% AnnotationDbi::columns(getpckg))
			stop(x, "is not a valid ID,check by loading library(" ,package,") then, AnnotationDbi::columns(",package,")")
	})

	return(NULL)
}



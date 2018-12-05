#' Get download URL
#'
#' Convert an output file name and location to a URL that can be used to
#' download the file.
#'
#' @param file name of the output file
#' @param folder name of the directory in the output directory containing the
#' output file
#'
#' @return URL linking to the file
getDownloadURL <- function(file, folder = NULL) {
    remote <- workflowr::wflow_git_remote(verbose = FALSE)["origin"]

    url <- gsub(":", "/", remote)
    url <- gsub("git@", "http://", url)
    url <- gsub(".git", "", url, fixed = TRUE)
    url <- paste(url, "raw/master/output", sep = "/")

    if (is.null(folder)) {
        url <- paste(url, file, sep = "/")
    } else {
        url <- paste(url, folder, file, sep = "/")
    }

    return(url)
}

#' Write gene table
#'
#' Save a gene table with the results of differential expression analysis to a
#' file
#'
#' @param gene.table gene table to save
#' @param path file path to save location
#'
#' @details
#' data.frame objects will be saved as a (zipped) CSV file. List objects will
#' be saved in XLSX format.
#'
#' @return URL linking to the file
writeGeneTable <- function(gene.table, path) {
    annot <- readr::read_tsv(here("data/biomart_annotation.tsv"),
                             col_types = readr::cols(
                                 .default = readr::col_character(),
                                 gc_content = readr::col_double()
                             ))
    annot <- dplyr::select(annot, feature_symbol, ensembl_id, entrezgene)

    if (is.data.frame(gene.table)) {
        gene.table <- dplyr::left_join(gene.table, annot,
                                       by = c(gene = "feature_symbol"))
        gene.table <- dplyr::select(gene.table, gene, ensembl_id, entrezgene,
                                    dplyr::everything())
        zip.path <- paste0(path, ".zip")
        if (file.exists(zip.path)) {file.remove(zip.path)}
        readr::write_csv(gene.table, path, na = "")
        zip(zip.path, path, flags = "-q -j")
        invisible(file.remove(path))
    } else {
        gene.table <- lapply(gene.table, function(x) {
            x <- dplyr::left_join(x, annot, by = c(Gene = "feature_symbol"))
            x <- dplyr::select(x, Gene, ENSEMBL = ensembl_id,
                               EntrezGene = entrezgene,
                               dplyr::everything())
        })
        writexl::write_xlsx(gene.table, path)
    }
}

#' Load TSV as SingleCellExperiment
#'
#' Loads 10x scRNA-seq data available as a TSV into a SingleCellExperiment
#' object
#'
#' @param dir Directory containing Cell Ranger counts matrix
#' @param genes data.frame where first column is gene ids and second column is
#'        gene symbols
#' @param dataset Name of the dataset
#' @param add.anno Whether to add feature annotation?
#' @param calc.qc Whether to calculate scater QC metrics?
#' @param calc.cpm Whether to calculate Counts Per Million?
#' @param pct.mt Whether to calculate percent mitochondrial?
#' @param pct.ribo Whether to calculate percent ribosomal?
#' @param cell.cycle Whether to assign cell cycle?
#' @param sparse Whethere to store data as a sparse matrix?
#' @param bpparam A BiocParallelParam object to use for parallel processing
#' @param verbose Whether to print progress messages?
#'
#' @return SingleCellExperiment object containing 10x data
loadTSVSCE <- function(path, genes, dataset, add.anno = FALSE, calc.qc = FALSE,
                       calc.cpm = FALSE, pct.mt = FALSE, pct.ribo = FALSE,
                       cell.cycle = FALSE, sparse = TRUE,
                       bpparam = BiocParallel::SerialParam(),
                       verbose = FALSE) {

    if (verbose) {message("Reading matrix...")}
    mtx <- as.matrix(read.delim(path, row.names = 1))

    cell.names <- paste0("Cell", 1:ncol(mtx))
    colnames(mtx) <- cell.names

    gene.vec <- genes[, 1]
    gene.vec.names <- make.unique(genes[, 2])
    names(gene.vec) <- gene.vec.names

    row.data <- data.frame(gene_id = gene.vec[rownames(mtx)],
                           symbol = names(gene.vec[rownames(mtx)]),
                           stringsAsFactors = FALSE)
    rownames(row.data) <- row.data$gene_id

    rownames(mtx) <- NULL

    if (sparse) {
        mtx <- as(mtx, "dgCMatrix")
    }

    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mtx),
                                                      rowData = row.data)

    SummarizedExperiment::colData(sce)$Cell <- cell.names
    SummarizedExperiment::colData(sce)$Dataset <- dataset
    colnames(SummarizedExperiment::rowData(sce)) <- c("gene_id", "symbol")

    sce <- annotateSCE(sce, org, add.anno = add.anno, calc.qc = calc.qc,
                       calc.cpm = calc.cpm, pct.mt = pct.mt,
                       pct.ribo = pct.ribo, cell.cycle = cell.cycle,
                       bpparam = bpparam, verbose = verbose)

    if (verbose) {message("Done!")}
    return(sce)
}


#' Load 10x data
#'
#' Loads 10x data into a SingleCellExperiment object
#'
#' @param dir Directory containing Cell Ranger counts matrix
#' @param dataset Name of the dataset
#' @param org Whether it is human or mouse data
#' @param add.anno Whether to add feature annotation?
#' @param calc.qc Whether to calculate QC metrics?
#' @param calc.cpm Whether to calculate Counts Per Million?
#' @param pct.mt Whether to calculate percent mitochondrial?
#' @param pct.ribo Whether to calculate percent ribosomal?
#' @param cell.cycle Whether to assign cell cycle?
#' @param kid.genes Whether to add kidney marker genes?
#' @param kid.path Path to kidney marker gene TSV file
#' @param kid.col Column containing ENSEMBL ids
#' @param cr.clusts Whether to add Cell Ranger clusters?
#' @param cr.tsne Whether to add Cell Ranger tSNE projection?
#' @param cr.genes Whether to add Cell Ranger DE results?
#' @param sparse Whethere to store data as a sparse matrix?
#' @param bpparam A BiocParallelParam object to use for parallel processing
#' @param verbose Whether to print progress messages?
#'
#' @return SingleCellExperiment object containing 10x data
load10xSCE <- function(dir, dataset, org = c("human", "mouse"),
                       add.anno = FALSE, calc.qc = FALSE, calc.cpm = FALSE,
                       pct.mt = FALSE, pct.ribo = FALSE, cell.cycle = FALSE,
                       sparse = TRUE, bpparam = BiocParallel::SerialParam(),
                       verbose = FALSE) {

    if (verbose) {message("Reading 10x data...")}
    sce <- DropletUtils::read10xCounts(dir, col.names = TRUE)

    if (!sparse) {
        BiocGenerics::counts(sce) <- as.matrix(BiocGenerics::counts(sce))
    }

    cell.names <- paste0("Cell", 1:ncol(sce))
    colnames(sce) <- cell.names

    barcodes <- SummarizedExperiment::colData(sce)$Barcode
    SummarizedExperiment::colData(sce)$Cell <- cell.names
    SummarizedExperiment::colData(sce)$Barcode <- stringr::str_sub(barcodes,
                                                                   end = -3)
    SummarizedExperiment::colData(sce)$Sample <- stringr::str_sub(barcodes,
                                                                  start = -1)
    SummarizedExperiment::colData(sce)$Dataset <- dataset
    colnames(SummarizedExperiment::rowData(sce)) <- c("gene_id", "symbol")

    sce <- annotateSCE(sce, org, add.anno = add.anno, calc.qc = calc.qc,
                       calc.cpm = calc.cpm, pct.mt = pct.mt,
                       pct.ribo = pct.ribo, cell.cycle = cell.cycle,
                       bpparam = bpparam, verbose = verbose)

    if (verbose) {message("Done!")}

    return(sce)
}


#' Annotate SCE
#'
#' Add additional metadata to a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object
#' @param org Whether it is human or mouse data
#' @param dataset Name of the dataset
#' @param add.anno Whether to add feature annotation?
#' @param calc.qc Whether to calculate QC metrics?
#' @param calc.cpm Whether to calculate Counts Per Million?
#' @param pct.mt Whether to calculate percent mitochondrial?
#' @param pct.ribo Whether to calculate percent ribosomal?
#' @param cell.cycle Whether to assign cell cycle?
#' @param bpparam A BiocParallelParam object to use for parallel processing
#' @param verbose Whether to print progress messages?
#'
#' @return SingleCellExperiment object with additional annotation
annotateSCE <- function(sce, org = c("human", "mouse"), add.anno = FALSE,
                        calc.qc = FALSE, calc.cpm = FALSE, pct.mt = FALSE,
                        pct.ribo = FALSE, cell.cycle = FALSE,
                        bpparam = BiocParallel::SerialParam(),
                        verbose = FALSE) {

    if (add.anno) {
        if (verbose) {message("Adding feature annotation...")}
        sce <- addFeatureAnnos(sce, org)
    }

    if (cell.cycle) {
        if (verbose) {message("Assigning cell cycle phases...")}
        sce <- addCellCycle(sce, org, bpparam, verbose)
    }

    if (calc.cpm) {
        if (verbose) {message("Calculating CPM...")}
        SingleCellExperiment::cpm(sce) <- scater::calculateCPM(sce,
                                                               use_size_factors = FALSE)
    }

    if (calc.qc) {
        if (verbose) {message("Calculating QC metrics...")}
        sce <- scater::calculateQCMetrics(sce)
        SummarizedExperiment::colData(sce)$pct_dropout <- 100 *
            (1 - SummarizedExperiment::colData(sce)$total_features / nrow(sce))
    }

    if (pct.mt) {
        if (verbose) {"Adding percentage mitochondrial..."}
        sce <- addPctCountsMT(sce)
    }

    if (pct.ribo) {
        if (verbose) {"Adding percentage ribosomal..."}
        sce <- addPctCountsRibo(sce)
    }

    return(sce)
}


#' Add feature annotations
#'
#' Add feature annotations to an SCE object. Just calls the
#' \code{getBMFeatureAnnos} function in \code{scater} with a specific set of
#' parameters.
#'
#' @param sce SCE object to add annotations to.
#' @param org Whether it is human or mouse data
#'
#' @return SCE with addtional feature annotations.
addFeatureAnnos <- function(sce, org) {

    if (org == "human") {
        dataset <- "hsapiens_gene_ensembl"
        symbol <- "hgnc_symbol"
    } else if (org == "mouse") {
        dataset <- "mmusculus_gene_ensembl"
        symbol <- "mgi_symbol"
    } else {
        stop("org should be 'human' or 'mouse'")
    }

    sce <- scater::getBMFeatureAnnos(sce,
                                     filters = "ensembl_gene_id",
                                     attributes = c("ensembl_gene_id",
                                                    "entrezgene",
                                                    "external_gene_name",
                                                    symbol,
                                                    "chromosome_name",
                                                    "description",
                                                    "gene_biotype",
                                                    "percentage_gene_gc_content"),
                                     feature_symbol = "external_gene_name",
                                     dataset = dataset)

    SummarizedExperiment::rowData(sce)$description <- gsub("\\s\\[.*\\]", "",
                                                           SummarizedExperiment::rowData(sce)$description)

    return(sce)
}


#' Add sample metadata
#'
#' Add sample metadata to an SCE object.
#'
#' @param sce SCE object to add metadata to.
#' @param dir Directory where processed data is stored.
#' @param cols Columns of the metadata file to add. If \code{NULL} than all
#'        columns with be added.
#'
#' @return SCE with additional sample metadata.
addSampleMeta <- function(sce, dir,
                          cols = c("BioSpecimenID", "Concentration", "TotalDNA",
                                   "DilutionConc", "LibraryPrepBatch")) {

    targets <- suppressWarnings(
        readr::read_tsv(file.path(dir, "metadata/SampleInfo.txt"),
                        col_types = readr::cols(
                            .default          = readr::col_character(),
                            Sample            = readr::col_integer(),
                            OriginalPlate     = readr::col_integer(),
                            Concentration     = readr::col_double(),
                            ApproxVolume      = readr::col_integer(),
                            TotalDNA          = readr::col_double(),
                            VolSampleRequired = readr::col_integer(),
                            LibraryPrepBatch  = readr::col_integer()
                        ))
    )

    if (!is.null(cols)) {
        if (!("BioSpecimenID" %in% cols)) {
            cols <- c("BioSpecimenID", cols)
        }
        targets <- dplyr::select_(targets, .dots = cols)
    }

    pheno <- dplyr::as_tibble(SummarizedExperiment::colData(sce))
    pheno <- dplyr::mutate(pheno,
                           BioSpecimenID = paste(Organoid, Well, sep = "_"))
    pheno <- dplyr::left_join(pheno, targets, by = "BioSpecimenID")
    rownames(pheno) <- pheno$Cell

    SummarizedExperiment::colData(sce) <- pheno

    return(sce)
}


#' Add Percent Counts Mitochondrial
#'
#' Add the percentage of counts in each cell in an SCE object that are
#' assigned to mitochondrial genes.
#'
#' @param sce DESCRIPTION.
#'
#' @details
#' Assumes \code{colData} has columns called \code{chomosome_name} and
#' \code{description} such as those created by \code{addFeatureAnnos}.
#'
#' Mitochondrial genes are defined as those on the MT chromosome or with
#' "mitochondrial" in the description.
#'
#' @return SCE with addtional "PctCountsMT" colData column
addPctCountsMT <- function(sce) {
    is.mt <- SummarizedExperiment::rowData(sce)$chromosome_name == "MT"
    is.mt[is.na(is.mt)] <- FALSE
    is.mt <- is.mt | grepl("mitochondrial", SummarizedExperiment::rowData(sce)$description)

    SummarizedExperiment::colData(sce)$PctCountsMT <- colSums(
        as.matrix(BiocGenerics::counts(sce)[is.mt, ])) /
        SummarizedExperiment::colData(sce)$total_counts * 100

    return(sce)
}

#' Add Percent Counts Ribosomal
#'
#' Add the percentage of counts in each cell in an SCE object that are
#' assigned to ribosomal genes.
#'
#' @param sce DESCRIPTION.
#'
#' @details
#' Assumes \code{colData} has a column called \code{description} such as that
#' created by \code{addFeatureAnnos}.
#'
#' Ribosomal genes are defined as those with "ribosom" in the description.
#'
#' @return SCE with addtional "PctCountsRibo" colData column
addPctCountsRibo <- function(sce) {
    is.ribo <- grepl("ribosom", SummarizedExperiment::rowData(sce)$description)

    SummarizedExperiment::colData(sce)$PctCountsRibo <- colSums(
        as.matrix(BiocGenerics::counts(sce)[is.ribo, ])) /
        SummarizedExperiment::colData(sce)$total_counts * 100

    return(sce)
}

#' Add cell cycle
#'
#' Add cell cycle phase to cells in an SCE object. Calculated using the
#' \code{\link[scran]{cyclone}} function.
#'
#' @param sce SCE object.
#' @param org Whether it is human or mouse data
#' @param bpparam A BiocParallelParam object to use for parallel processing
#' @param verbose Whether to print progress messages?
#'
#' @return SCE object with assigned cell cycles
addCellCycle <- function(sce, org, bpparam, verbose) {

    if (org == "human") {
        cc.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds",
                                        package = "scran"))
    } else {
        cc.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds",
                                        package = "scran"))
    }

    cycles <- scran::cyclone(sce, pairs = cc.pairs, BPPARAM = bpparam,
                             verbose = verbose)
    SummarizedExperiment::colData(sce)$G1Score <- cycles$scores$G1
    SummarizedExperiment::colData(sce)$SScore <- cycles$scores$G1
    SummarizedExperiment::colData(sce)$G2MScore <- cycles$scores$G2M
    SummarizedExperiment::colData(sce)$CellCycle <- cycles$phases

    return(sce)
}

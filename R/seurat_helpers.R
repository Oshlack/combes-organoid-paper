#' Prefix cell names
#'
#' Add a prefix to cell names in a Seurat object. Useful before merging.
#'
#' @param object Seurat object to add prefix to
#' @param prefix prefix to add
#'
#' @return Seurat object with new cell names
prefixCellNames <- function(object, prefix) {

    object@cell.names <- paste(prefix ,object@cell.names, sep = "_")
    colnames(object@raw.data) <- paste(
        prefix, colnames(object@raw.data), sep = "_"
    )
    rownames(object@meta.data) <- paste(
        prefix, rownames(object@meta.data), sep = "_"
    )

    return(object)
}

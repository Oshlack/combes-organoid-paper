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
    remote <- workflowr::wflow_git_remote(verbose = FALSE)

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

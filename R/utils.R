#' @importFrom utils download.file head packageVersion
#' @importFrom stats setNames
NULL

`%||%` <- function(a, b) if (is.null(a)) b else a

## TRUE if `x` is missing, empty or only contains blank/NA values.
.is_blank <- function(x) {
    is.null(x) || !length(x) || all(is.na(x) | !nzchar(trimws(x)))
}

## First non-blank value of a (possibly NULL) element, or `default`.
.first <- function(x, default = NA_character_) {
    if (.is_blank(x)) default else trimws(as.character(x[[1L]]))
}

## Remove characters that are not allowed within a cell of the tab-delimited
## mzTab-M file.
.clean_text <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    trimws(gsub("[\t\r\n]+", " ", x))
}

## Quote a name or value used within a CV parameter if it contains characters
## that would break the `[label, accession, name, value]` syntax.
.cv_quote <- function(x) {
    x <- gsub("\"", "'", .clean_text(x), fixed = TRUE)
    quote <- grepl("[],[]", x)
    x[quote] <- paste0("\"", x[quote], "\"")
    x
}

## A user parameter, i.e. a CV parameter without CV label and accession.
.user_param <- function(name, value = "") {
    paste0("[,, ", .cv_quote(name), ", ", .cv_quote(value), "]")
}

.cv_param <- function(label, accession, name, value = "") {
    paste0("[", label, ", ", accession, ", ", name, ", ", .cv_quote(value), "]")
}

.pkg_version <- function() {
    as.character(packageVersion("mwtab2mztabm"))
}

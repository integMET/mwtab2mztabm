.MW_REST <- "https://www.metabolomicsworkbench.org/rest/study"

#' Locate and download mwTab JSON files from the Metabolomics Workbench
#'
#' `mwtabUrl()` resolves a Metabolomics Workbench study ID (e.g.
#' `"ST003934"`), analysis ID (e.g. `"AN006459"`) or study web page URL to the
#' URL of the REST service returning the corresponding mwTab JSON file. URLs
#' of study pages (such as
#' `https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST003934`)
#' return HTML; for these, the study or analysis ID is extracted from the URL
#' and the REST URL is returned instead. Any other URL is returned as is.
#'
#' `downloadMwTab()` downloads the mwTab JSON file and checks that the
#' response is indeed JSON.
#'
#' A study-level request returns all analyses of the study within a single
#' JSON file.
#'
#' @param x `character(1)` with a study ID, analysis ID or URL.
#'
#' @param destdir `character(1)` with the directory to save the file to.
#'
#' @param overwrite `logical(1)` whether an existing file should be replaced.
#'
#' @param quiet `logical(1)` passed to [utils::download.file()].
#'
#' @return `mwtabUrl()`: `character(1)` with the URL. `downloadMwTab()`:
#'     `character(1)` with the path to the downloaded file.
#'
#' @export
#'
#' @examples
#' mwtabUrl("ST003934")
#' mwtabUrl("AN006459")
#' mwtabUrl(paste0("https://www.metabolomicsworkbench.org/data/",
#'     "DRCCMetadata.php?Mode=Study&StudyID=ST003934&StudyType=MS&ResultType=1"))
#'
#' \donttest{
#' f <- downloadMwTab("ST003934")
#' }
mwtabUrl <- function(x) {
    if (!is.character(x) || length(x) != 1L || is.na(x))
        stop("'x' has to be a single study ID, analysis ID or URL",
             call. = FALSE)
    x <- trimws(x)
    if (grepl("^ST[0-9]+$", x, ignore.case = TRUE))
        return(paste0(.MW_REST, "/study_id/", toupper(x), "/mwtab/json"))
    if (grepl("^AN[0-9]+$", x, ignore.case = TRUE))
        return(paste0(.MW_REST, "/analysis_id/", toupper(x), "/mwtab/json"))
    if (!grepl("^(https?|ftp|file)://", x))
        stop("'", x, "' is neither a Metabolomics Workbench study/analysis ",
             "ID nor a URL", call. = FALSE)
    if (grepl("metabolomicsworkbench\\.org", x) && !grepl("/rest/", x)) {
        an <- .url_id(x, "ANALYSIS_?ID", "AN")
        if (!is.na(an)) return(mwtabUrl(an))
        st <- .url_id(x, "STUDY_?ID", "ST")
        if (!is.na(st)) return(mwtabUrl(st))
    }
    x
}

## Extract a `<prefix>000000` ID from the query parameter `param` of `url`.
.url_id <- function(url, param, prefix) {
    m <- regmatches(url, regexec(
        paste0("[?&]", param, "=(", prefix, "[0-9]+)"), url,
        ignore.case = TRUE))[[1L]]
    if (length(m)) toupper(m[2L]) else NA_character_
}

#' @rdname mwtabUrl
#'
#' @export
downloadMwTab <- function(x, destdir = tempdir(), overwrite = FALSE,
                          quiet = TRUE) {
    url <- mwtabUrl(x)
    id <- regmatches(url, regexpr("(ST|AN)[0-9]+", url))
    if (!length(id)) id <- "mwtab"
    dest <- file.path(destdir, paste0(id, "_mwtab.json"))
    if (file.exists(dest) && !overwrite)
        return(dest)
    dir.create(destdir, recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(fileext = ".json")
    on.exit(unlink(tmp))
    download.file(url, tmp, mode = "wb", quiet = quiet)
    if (!.is_json_file(tmp))
        stop("The URL '", url, "' did not return a JSON file. Please provide ",
             "a Metabolomics Workbench study ID, analysis ID or the URL of ",
             "an mwTab JSON file.", call. = FALSE)
    if (!file.copy(tmp, dest, overwrite = TRUE))
        stop("Failed to save the mwTab file to '", dest, "'", call. = FALSE)
    dest
}

.is_json_file <- function(f) {
    start <- readChar(f, 256L, useBytes = TRUE)
    length(start) && grepl("^[[:space:]]*[[{]", start, useBytes = TRUE)
}

#' Read an mwTab JSON file
#'
#' `readMwTab()` reads an mwTab JSON file (as returned by the Metabolomics
#' Workbench REST service) and returns its analyses. Each analysis is a named
#' `list` with the mwTab sections (`"METABOLOMICS WORKBENCH"`, `"PROJECT"`,
#' `"STUDY"`, `"SUBJECT_SAMPLE_FACTORS"`, `"MS_METABOLITE_DATA"`, ...).
#'
#' @param file `character(1)` with the path to the mwTab JSON file.
#'
#' @return Named `list` of analyses, names being the analysis IDs.
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "example_mwtab.json", package = "mwtab2mztabm")
#' mw <- readMwTab(f)
#' names(mw)
#' names(mw[[1L]])
readMwTab <- function(file) {
    x <- jsonlite::fromJSON(file, simplifyVector = FALSE)
    if (.is_mwtab_analysis(x))
        x <- list(x)
    if (!length(x) || !all(vapply(x, .is_mwtab_analysis, NA)))
        stop("'", file, "' is not an mwTab JSON file", call. = FALSE)
    ids <- vapply(x, function(z)
        .first(z[["METABOLOMICS WORKBENCH"]][["ANALYSIS_ID"]]), "")
    ids[is.na(ids)] <- paste0("analysis_", which(is.na(ids)))
    names(x) <- ids
    x
}

.is_mwtab_analysis <- function(x) {
    is.list(x) && !is.null(names(x)) &&
        is.list(x[["METABOLOMICS WORKBENCH"]])
}

## Return a single mwTab analysis from `x`, which can be the path to an mwTab
## JSON file, the result of `readMwTab()` or a single analysis.
.mwtab_analysis <- function(x, analysis_id = NULL) {
    if (is.character(x) && length(x) == 1L)
        x <- readMwTab(x)
    if (.is_mwtab_analysis(x))
        return(x)
    if (!is.list(x) || !length(x))
        stop("'x' has to be an mwTab JSON file or the result of readMwTab()",
             call. = FALSE)
    if (is.null(analysis_id)) {
        if (length(x) > 1L)
            stop("The mwTab file contains ", length(x), " analyses (",
                 paste(names(x), collapse = ", "), "). Please select one ",
                 "with parameter 'analysis_id'.", call. = FALSE)
        return(x[[1L]])
    }
    if (!analysis_id %in% names(x))
        stop("Analysis '", analysis_id, "' not found. Available analyses: ",
             paste(names(x), collapse = ", "), call. = FALSE)
    x[[analysis_id]]
}

## Name of the section with the metabolite data of an analysis.
.data_section_names <- function(x) {
    grep("_(METABOLITE|BINNED)_DATA$", names(x), value = TRUE)
}

## Sample factors of an analysis as a list of entries with elements
## `"Subject ID"`, `"Sample ID"`, `"Factors"` and `"Additional sample data"`,
## the latter two being named lists of `character(1)`.
.mwtab_sample_factors <- function(x) {
    lapply(x[["SUBJECT_SAMPLE_FACTORS"]], function(z) {
        list(`Subject ID` = .first(z[["Subject ID"]], ""),
             `Sample ID` = .first(z[["Sample ID"]], ""),
             Factors = .key_values(z[["Factors"]]),
             `Additional sample data` =
                 .key_values(z[["Additional sample data"]]))
    })
}

## Factors are a JSON object in current mwTab JSON files, but might also be
## provided as `"key1:value1 | key2:value2"` string.
.key_values <- function(x) {
    if (is.null(x))
        return(setNames(list(), character()))
    if (is.character(x)) {
        parts <- trimws(strsplit(paste(x, collapse = "|"), "|",
                                 fixed = TRUE)[[1L]])
        parts <- parts[nzchar(parts)]
        keys <- trimws(sub(":.*$", "", parts))
        vals <- trimws(ifelse(grepl(":", parts), sub("^[^:]*:", "", parts), ""))
        x <- setNames(as.list(vals), keys)
    }
    lapply(x, function(v) if (is.null(v)) "" else paste(unlist(v),
                                                         collapse = "; "))
}

#' Export and import the mwTab sample metadata as TOML
#'
#' `mwtabToToml()` writes the metadata of one mwTab analysis to a TOML file.
#' All mwTab sections are exported with their original names and keys
#' (`"METABOLOMICS WORKBENCH"`, `"PROJECT"`, `"STUDY"`, `"SUBJECT"`,
#' `"COLLECTION"`, `"MS"`, ...) except for the metabolite data itself, of which
#' only the `Units` are kept. The sample metadata (mwTab section
#' `SUBJECT_SAMPLE_FACTORS`) is exported as an array of tables, one per sample,
#' with sub-tables `Factors` and `"Additional sample data"`.
#'
#' Since TOML tables are unordered, the original order of the sample factors
#' and the additional sample data is stored in the table `[mwtab2mztabm]`.
#'
#' `readMwTabToml()` reads such a TOML file.
#'
#' @param x for `mwtabToToml()`: the path to an mwTab JSON file, the result of
#'     [readMwTab()] or a single mwTab analysis.
#'
#' @param file `character(1)` with the path of the TOML file.
#'
#' @param analysis_id `character(1)` with the ID of the analysis to export.
#'     Only needed if the mwTab file contains more than one analysis.
#'
#' @return `mwtabToToml()`: `file`, invisibly. `readMwTabToml()`: a named
#'     `list` with the content of the TOML file.
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "example_mwtab.json", package = "mwtab2mztabm")
#' toml <- tempfile(fileext = ".toml")
#' mwtabToToml(f, toml, analysis_id = "AN999998")
#' cat(head(readLines(toml), 30), sep = "\n")
#'
#' meta <- readMwTabToml(toml)
#' meta$SUBJECT_SAMPLE_FACTORS[[1L]]
mwtabToToml <- function(x, file, analysis_id = NULL) {
    x <- .mwtab_analysis(x, analysis_id)
    data_sections <- .data_section_names(x)
    meta <- x[!names(x) %in% c(data_sections, "SUBJECT_SAMPLE_FACTORS")]
    for (s in data_sections)
        meta[[s]] <- list(Units = .first(x[[s]][["Units"]], ""))
    ssf <- .mwtab_sample_factors(x)
    ord <- list(
        format_version = 1L,
        factor_order = .key_order(ssf, "Factors"),
        additional_data_order = .key_order(ssf, "Additional sample data"))
    mw <- x[["METABOLOMICS WORKBENCH"]]
    header <- c(
        paste0("# mwTab metadata of study ", .first(mw$STUDY_ID, "?"),
               ", analysis ", .first(mw$ANALYSIS_ID, "?")),
        paste0("# Exported by mwtab2mztabm ", .pkg_version(), "."),
        "# Sample metadata is listed in [[SUBJECT_SAMPLE_FACTORS]].")
    lines <- c(header, "",
               .toml_table(c(list(mwtab2mztabm = ord), meta,
                             list(SUBJECT_SAMPLE_FACTORS = ssf))))
    con <- file(file, open = "wb")
    on.exit(close(con))
    writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
    invisible(file)
}

.key_order <- function(ssf, element) {
    unique(unlist(lapply(ssf, function(z) names(z[[element]])),
                  use.names = FALSE)) %||% character()
}

#' @rdname mwtabToToml
#'
#' @export
readMwTabToml <- function(file) {
    x <- unclass(RcppTOML::parseTOML(file, escape = FALSE))
    attr(x, "file") <- NULL
    ord <- x[["mwtab2mztabm"]]
    x[["SUBJECT_SAMPLE_FACTORS"]] <- lapply(
        x[["SUBJECT_SAMPLE_FACTORS"]], function(z) {
            z[["Factors"]] <- .reorder(z[["Factors"]], ord$factor_order)
            z[["Additional sample data"]] <- .reorder(
                z[["Additional sample data"]], ord$additional_data_order)
            z
        })
    x
}

.reorder <- function(x, ord) {
    if (!length(x)) return(list())
    x[c(intersect(ord, names(x)), setdiff(names(x), ord))]
}

## Serialization of R lists to TOML. Named lists are written as tables,
## unnamed lists of named lists as arrays of tables and everything else as
## key/value pairs.
.toml_table <- function(x, path = character()) {
    x <- x[!vapply(x, is.null, NA)]
    is_tbl <- vapply(x, .is_toml_table, NA)
    is_aot <- vapply(x, .is_toml_array_of_tables, NA)
    keys <- names(x)
    out <- character()
    for (k in keys[!is_tbl & !is_aot])
        out <- c(out, paste0(.toml_key(k), " = ", .toml_value(x[[k]])))
    for (k in keys[is_tbl]) {
        p <- c(path, k)
        out <- c(out, if (length(out)) "", paste0("[", .toml_path(p), "]"),
                 .toml_table(x[[k]], p))
    }
    for (k in keys[is_aot]) {
        p <- c(path, k)
        for (el in x[[k]])
            out <- c(out, if (length(out)) "", paste0("[[", .toml_path(p), "]]"),
                     .toml_table(el, p))
    }
    out
}

.is_toml_table <- function(x) {
    is.list(x) && !is.null(names(x))
}

.is_toml_array_of_tables <- function(x) {
    is.list(x) && is.null(names(x)) && length(x) &&
        all(vapply(x, .is_toml_table, NA))
}

.toml_key <- function(x) {
    bare <- grepl("^[A-Za-z0-9_-]+$", x)
    x[!bare] <- .toml_string(x[!bare])
    x
}

.toml_path <- function(x) paste(.toml_key(x), collapse = ".")

.toml_value <- function(x) {
    if (is.list(x))
        x <- unlist(x, use.names = FALSE)
    if (!length(x))
        return("[]")
    vals <- if (is.logical(x)) {
        ifelse(is.na(x), "\"\"", tolower(as.character(x)))
    } else if (is.integer(x)) {
        ifelse(is.na(x), "\"\"", as.character(x))
    } else if (is.numeric(x)) {
        ifelse(is.na(x), "nan", ifelse(is.infinite(x),
                                       ifelse(x > 0, "inf", "-inf"),
                                       format(x, digits = 15L, trim = TRUE)))
    } else .toml_string(x)
    if (length(vals) == 1L) vals
    else paste0("[", paste(vals, collapse = ", "), "]")
}

.TOML_ESCAPES <- c("\\" = "\\\\", "\"" = "\\\"", "\b" = "\\b", "\t" = "\\t",
                   "\n" = "\\n", "\f" = "\\f", "\r" = "\\r")

.toml_string <- function(x) {
    x <- enc2utf8(as.character(x))
    x[is.na(x)] <- ""
    m <- gregexpr("[\\\\\"\001-\037\177]", x, perl = TRUE)
    regmatches(x, m) <- lapply(regmatches(x, m), function(ch) {
        esc <- unname(.TOML_ESCAPES[ch])
        miss <- is.na(esc)
        esc[miss] <- sprintf("\\u%04X", vapply(ch[miss], utf8ToInt, 1L))
        esc
    })
    paste0("\"", x, "\"")
}

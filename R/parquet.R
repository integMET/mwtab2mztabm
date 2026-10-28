#' Export and import the mwTab metabolite quantification table as Parquet
#'
#' `mwtabToParquet()` writes the metabolite-wise quantification table of one
#' mwTab analysis (section `MS_METABOLITE_DATA`) to a Parquet file. Each row is
#' one metabolite. The columns are:
#'
#' - `Metabolite`: the metabolite name as reported in mwTab.
#' - the metabolite annotations from the `Metabolites` block of the mwTab
#'   section (e.g. `"pubchem_id"`, `"kegg_id"`, `"HMDB"`, `"retention index"`),
#'   matched to the rows of the quantification table by metabolite name.
#'   Names without an exact match are matched after normalizing the
#'   characters the Metabolomics Workbench replaces in the quantification
#'   table (`"^"` for `"'"`, `"_"` for `","` and the numbering of repeated
#'   names), if that match is unique.
#' - one `double` column with the abundances for each sample, in the order of
#'   the mwTab quantification table. Values that can not be converted to a
#'   number (such as `""` or `"<LOD"`) are stored as missing values.
#'
#' The names of the sample and annotation columns, the units of the
#' abundances and the study and analysis IDs are stored as key-value metadata
#' of the Parquet file (keys `"mwtab2mztabm.sample_columns"`,
#' `"mwtab2mztabm.annotation_columns"`, `"mwtab.units"`, `"mwtab.study_id"`
#' and `"mwtab.analysis_id"`).
#'
#' `readMwTabParquet()` reads such a Parquet file.
#'
#' @inheritParams mwtabToToml
#'
#' @param file `character(1)` with the path of the Parquet file.
#'
#' @param sample_columns optional `character` with the names of the columns
#'     containing the abundances. If not provided, they are taken from the
#'     metadata of the Parquet file or, if not available there, all numeric
#'     columns are used.
#'
#' @return `mwtabToParquet()`: `file`, invisibly. `readMwTabParquet()`: a
#'     `data.frame` with the content of the Parquet file and attributes
#'     `"sample_columns"`, `"annotation_columns"` and `"metadata"` (named
#'     `character` with the key-value metadata of the file).
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "example_mwtab.json", package = "mwtab2mztabm")
#' pq <- tempfile(fileext = ".parquet")
#' mwtabToParquet(f, pq, analysis_id = "AN999998")
#'
#' tab <- readMwTabParquet(pq)
#' tab
#' attr(tab, "sample_columns")
mwtabToParquet <- function(x, file, analysis_id = NULL) {
    x <- .mwtab_analysis(x, analysis_id)
    sec <- .data_section_names(x)
    if (!length(sec))
        stop("The mwTab analysis does not contain metabolite data (no ",
             "'MS_METABOLITE_DATA' section). The results might be ",
             "provided in a separate file.", call. = FALSE)
    d <- x[[sec[1L]]]
    dat <- .rows_to_df(d[["Data"]])
    if (!nrow(dat))
        stop("The '", sec[1L], "' section contains no data", call. = FALSE)
    name_col <- if ("Metabolite" %in% colnames(dat)) "Metabolite"
                else colnames(dat)[1L]
    samples <- setdiff(colnames(dat), name_col)
    abund <- lapply(dat[samples], .as_abundance)
    nbad <- sum(vapply(abund, attr, 1L, "n_invalid"))
    if (nbad)
        warning(nbad, " non-numeric abundance value(s) were replaced by ",
                "missing values", call. = FALSE)
    ann <- .metabolite_annotations(d[["Metabolites"]], dat[[name_col]])
    clash <- colnames(ann) %in% c("Metabolite", samples)
    colnames(ann)[clash] <- paste0("annotation_", colnames(ann)[clash])
    res <- data.frame(Metabolite = dat[[name_col]], ann,
                      lapply(abund, as.vector), check.names = FALSE)
    mw <- x[["METABOLOMICS WORKBENCH"]]
    meta <- c(
        mwtab2mztabm.sample_columns = .json(samples),
        mwtab2mztabm.annotation_columns = .json(colnames(ann)),
        mwtab.units = .first(d[["Units"]], ""),
        mwtab.study_id = .first(mw[["STUDY_ID"]], ""),
        mwtab.analysis_id = .first(mw[["ANALYSIS_ID"]], ""))
    nanoparquet::write_parquet(res, file, metadata = meta)
    invisible(file)
}

#' @rdname mwtabToParquet
#'
#' @export
readMwTabParquet <- function(file, sample_columns = NULL) {
    res <- as.data.frame(nanoparquet::read_parquet(file))
    kv <- nanoparquet::read_parquet_metadata(file)$file_meta_data$
        key_value_metadata[[1L]]
    meta <- setNames(kv$value, kv$key)
    meta <- meta[names(meta) != "ARROW:schema"]
    if (is.null(sample_columns)) {
        sample_columns <- if (is.na(meta["mwtab2mztabm.sample_columns"]))
            colnames(res)[vapply(res, is.numeric, NA)]
        else jsonlite::fromJSON(meta[["mwtab2mztabm.sample_columns"]])
    }
    if (!all(sample_columns %in% colnames(res)))
        stop("Sample column(s) ", paste0("\"", setdiff(sample_columns,
                                                      colnames(res)), "\"",
                                         collapse = ", "),
             " not found in '", file, "'", call. = FALSE)
    attr(res, "sample_columns") <- sample_columns
    attr(res, "annotation_columns") <-
        setdiff(colnames(res), c("Metabolite", sample_columns))
    attr(res, "metadata") <- meta
    res
}

.json <- function(x) {
    as.character(jsonlite::toJSON(as.character(x)))
}

## Convert a list of JSON objects (rows) into a character data.frame keeping
## the order of the keys.
.rows_to_df <- function(rows) {
    keys <- unique(unlist(lapply(rows, names), use.names = FALSE))
    cols <- lapply(keys, function(k) vapply(rows, function(r) {
        v <- r[[k]]
        if (is.null(v) || !length(v)) NA_character_ else as.character(v[[1L]])
    }, NA_character_))
    names(cols) <- keys
    as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
}

## Values considered as "missing" in mwTab quantification tables.
.NA_STRINGS <- c("", "NA", "N/A", "NaN", "nan", "null", "NULL", "-")

.as_abundance <- function(x) {
    x <- trimws(x)
    x[x %in% .NA_STRINGS] <- NA_character_
    res <- suppressWarnings(as.numeric(x))
    attr(res, "n_invalid") <- sum(is.na(res) & !is.na(x))
    res
}

## Annotation table aligned to the metabolites (`names`) of the
## quantification table.
.metabolite_annotations <- function(x, names) {
    empty <- data.frame(row.names = seq_along(names))
    if (!length(x))
        return(empty)
    ann <- .rows_to_df(x)
    key <- intersect(c("Metabolite", "metabolite_name"), colnames(ann))
    key <- if (length(key)) key[1L] else colnames(ann)[1L]
    ann_names <- ann[[key]]
    ann <- ann[, colnames(ann) != key, drop = FALSE]
    if (!ncol(ann))
        return(empty)
    if (length(ann_names) == length(names) &&
        identical(trimws(ann_names), trimws(names)))
        idx <- seq_along(names)
    else {
        idx <- match(trimws(names), trimws(ann_names))
        idx <- .match_normalized(names, ann_names, idx)
        if (anyNA(idx))
            warning(sum(is.na(idx)), " metabolite(s) have no annotation ",
                    "in the 'Metabolites' block", call. = FALSE)
    }
    ann <- ann[idx, , drop = FALSE]
    rownames(ann) <- NULL
    ann
}

## The Metabolomics Workbench writes some metabolite names differently in the
## quantification table and the `Metabolites` block: "'" as "^", "," as "_"
## and the suffixes of repeated names as "_1" instead of "_2". Metabolites
## without a match by name are matched to the still unused annotations by
## their normalized name, if that match is unique.
.match_normalized <- function(names, ann_names, idx) {
    miss <- which(is.na(idx))
    if (!length(miss))
        return(idx)
    norm <- function(x)
        sub("_[0-9]+$", "", chartr("^,", "'_", trimws(x)))
    free <- setdiff(seq_along(ann_names), idx)
    nfree <- norm(ann_names[free])
    nmiss <- norm(names[miss])
    ok <- nmiss %in% nfree[!duplicated(nfree) & !duplicated(nfree,
                                                             fromLast = TRUE)] &
        !nmiss %in% nmiss[duplicated(nmiss)]
    if (any(ok)) {
        idx[miss[ok]] <- free[match(nmiss[ok], nfree)]
        message(sum(ok), " metabolite(s) were matched to the 'Metabolites' ",
                "block by their normalized name")
    }
    idx
}

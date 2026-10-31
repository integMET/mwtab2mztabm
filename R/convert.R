#' Assemble an mzTab-M file from the TOML and Parquet files
#'
#' `assembleMzTabM()` creates the MTD section from the TOML file (with
#' [mtdFromToml()]) and the SML section from the Parquet file (with
#' [smlFromParquet()]), combines them into an [RmzTabM::MzTabM()] object and,
#' if `file` is provided, writes it with [RmzTabM::writeMzTabM()].
#'
#' The assays of the MTD section are defined in the order of the sample
#' columns of the Parquet file. Samples with metadata but without abundances
#' are not reported.
#'
#' The mzTab-M file is written in UTF-8 regardless of the locale of the R
#' session. If the session does not use a UTF-8 locale, the character type
#' (`LC_CTYPE`) is changed to a UTF-8 locale during this function. A warning is
#' issued if the system has no such locale, as non-ASCII characters might then be
#' corrupted.
#'
#' @param toml `character(1)` with the path to the TOML file created by
#'     [mwtabToToml()].
#'
#' @param parquet `character(1)` with the path to the Parquet file created by
#'     [mwtabToParquet()].
#'
#' @param file optional `character(1)` with the path of the mzTab-M file to
#'     write.
#'
#' @param validate `logical(1)` whether the written file should be validated
#'     with the jmzTab-M validator (see [RmzTabM::mzTabMValidator()]). The
#'     validator is downloaded on first use. The study variable abundances are
#'     also checked with [checkStudyVariableAbundance()]. An error is thrown if
#'     the validation fails; the invalid file is kept for inspection.
#'
#' @inheritParams mtdFromToml
#'
#' @inheritParams smlFromParquet
#'
#' @return The `MzTabM` object (invisibly if `file` was provided).
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "example_mwtab.json", package = "mwtab2mztabm")
#' toml <- tempfile(fileext = ".toml")
#' pq <- tempfile(fileext = ".parquet")
#' mwtabToToml(f, toml, analysis_id = "AN999998")
#' mwtabToParquet(f, pq, analysis_id = "AN999998")
#'
#' mzt <- assembleMzTabM(toml, pq)
#' mzt
#'
#' out <- tempfile(fileext = ".mztab")
#' assembleMzTabM(toml, pq, file = out, validate = FALSE)
#' cat(head(readLines(out), 10), sep = "\n")
assembleMzTabM <- function(toml, parquet, file = NULL, validate = TRUE,
                           reliability = "2", scan_polarity = NULL) {
    .with_utf8_ctype({
        tab <- readMwTabParquet(parquet)
        mtd <- mtdFromToml(toml, assays = attr(tab, "sample_columns"),
                           databases = smlDatabases(tab),
                           scan_polarity = scan_polarity)
        sml <- smlFromParquet(tab, mtd = mtd, reliability = reliability)
        mzt <- RmzTabM::MzTabM(mtd, sml = sml)
        if (!is.null(file)) {
            RmzTabM::writeMzTabM(mzt, path = file, validate = FALSE)
            if (validate) {
                .validate_mztab(file)
                .validate_study_variable_abundance(file)
            }
        }
    })
    if (is.null(file)) mzt else invisible(mzt)
}

## RmzTabM::mzTabMValidator() reports a valid file without any validation
## message as failed: it then counts lines containing "ERROR" in the output,
## which includes a logging configuration message of the validator. Such
## files are thus recognized here; all others are validated with RmzTabM.
.validate_mztab <- function(file) {
    out <- suppressWarnings(system2(
        RmzTabM::mzTabMValidator_download(), c("-c", shQuote(file)),
        stdout = TRUE, stderr = TRUE))
    if (any(grepl("There were 0 validation messages", out, fixed = TRUE)) &&
        any(grepl("No structural or logical errors", out, fixed = TRUE))) {
        message("Validation finished successfully: 0 message(s) found.")
        return(invisible(TRUE))
    }
    RmzTabM::mzTabMValidator(file)
    invisible(TRUE)
}

.validate_study_variable_abundance <- function(file) {
    res <- checkStudyVariableAbundance(file)
    if (nrow(res))
        stop("'", file, "': ", nrow(res), " abundance_study_variable value(s) ",
             "differ from the mean of the abundance_assay values. First ",
             "difference: SML_ID ", res$SML_ID[1L], ", ", res$study_variable[1L],
             " (reported ", res$reported[1L], ", expected ", res$expected[1L],
             ")", call. = FALSE)
    invisible(TRUE)
}

#' Check the study variable abundances of an mzTab-M file
#'
#' `checkStudyVariableAbundance()` reads an mzTab-M file and compares each
#' `abundance_study_variable[n]` value of the SML section with the mean of the
#' `abundance_assay[n]` values of the assays listed in
#' `study_variable[n]-assay_refs` of the MTD section.
#'
#' With the default `na.rm = TRUE` the mean of the available abundances is
#' expected, as calculated by [smlFromParquet()]: it is only missing (`null`) if
#' the abundances of all assays are missing. With `na.rm = FALSE` the mean is
#' expected to be missing if the abundance of any of the assays is missing.
#'
#' @param file `character(1)` with the path to the mzTab-M file.
#'
#' @param na.rm `logical(1)` whether missing assay abundances should be ignored
#'     when calculating the mean.
#'
#' @param tolerance `numeric(1)` with the accepted relative difference.
#'
#' @return `data.frame` with one row per difference and columns `"SML_ID"`,
#'     `"study_variable"`, `"reported"` and `"expected"`. It has no rows if all
#'     values are equal.
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "example_mwtab.json", package = "mwtab2mztabm")
#' res <- mwtabToMzTabM(f, outdir = tempdir(), analysis_id = "AN999998",
#'                      validate = FALSE)
#'
#' checkStudyVariableAbundance(res$mztab)
#' checkStudyVariableAbundance(res$mztab, na.rm = FALSE)
checkStudyVariableAbundance <- function(file, na.rm = TRUE, tolerance = 1e-8) {
    lines <- readLines(file)
    mtd <- .mztab_section(lines, "MTD")
    sml <- .mztab_section(lines, "SML")
    colnames(sml) <- .mztab_section(lines, "SMH")[1L, ]
    refs <- which(grepl("^study_variable\\[[0-9]+\\]-assay_refs$", mtd[, 2L]))
    res <- lapply(refs, function(i) {
        sv <- sub("-assay_refs$", "", mtd[i, 2L])
        svcol <- paste0("abundance_", sv)
        acols <- paste0("abundance_",
                        strsplit(mtd[i, 3L], "|", fixed = TRUE)[[1L]])
        miss <- setdiff(c(svcol, acols), colnames(sml))
        if (length(miss))
            stop("Columns ", paste(miss, collapse = ", "), " not found in the ",
                 "SML section of '", file, "'", call. = FALSE)
        reported <- .mztab_abundance(sml[, svcol])
        expected <- rowMeans(matrix(.mztab_abundance(sml[, acols]), nrow(sml)),
                             na.rm = na.rm)
        expected[is.nan(expected)] <- NA_real_
        ok <- is.na(reported) == is.na(expected)
        both <- !is.na(reported) & !is.na(expected)
        ok[both] <- abs(reported[both] - expected[both]) <=
            tolerance * abs(expected[both])
        data.frame(SML_ID = as.integer(sml[!ok, "SML_ID"]),
                   study_variable = rep(sv, sum(!ok)),
                   reported = reported[!ok], expected = expected[!ok])
    })
    res <- do.call(rbind, c(list(data.frame(
        SML_ID = integer(), study_variable = character(),
        reported = numeric(), expected = numeric())), res))
    rownames(res) <- NULL
    res
}

## Cells of the lines of a section of an mzTab-M file as a character matrix.
.mztab_section <- function(lines, prefix) {
    lines <- lines[startsWith(lines, paste0(prefix, "\t"))]
    if (!length(lines))
        stop("No ", prefix, " lines found in the mzTab-M file", call. = FALSE)
    do.call(rbind, strsplit(lines, "\t", fixed = TRUE))
}

.mztab_abundance <- function(x) {
    x[x == "null"] <- NA_character_
    as.numeric(x)
}

#' Convert mwTab files from the Metabolomics Workbench to mzTab-M
#'
#' `mwtabToMzTabM()` runs the full conversion for each analysis of an mwTab
#' file:
#'
#' 1. download the mwTab JSON file ([downloadMwTab()]), unless `x` is a local
#'    file,
#' 2. export the sample metadata to a TOML file ([mwtabToToml()]),
#' 3. export the metabolite quantification table to a Parquet file
#'    ([mwtabToParquet()]),
#' 4. create the MTD section from the TOML and the SML section from the
#'    Parquet file and write the mzTab-M file ([assembleMzTabM()]).
#'
#' The files are named `<STUDY_ID>_<ANALYSIS_ID>.toml`, `.parquet` and
#' `.mztab`. If the IDs can not be determined (see [readMwTab()]), the name of
#' the analysis (e.g. `analysis_1`) is used so that each analysis has its own
#' files.
#'
#' @param x `character(1)` with the path to a local mwTab JSON file or a
#'     Metabolomics Workbench study ID, analysis ID or URL (see
#'     [mwtabUrl()]).
#'
#' @param outdir `character(1)` with the directory for the output files.
#'
#' @param analysis_id optional `character` with the IDs of the analyses to
#'     convert. Defaults to all analyses in the mwTab file.
#'
#' @param ... additional parameters passed to [assembleMzTabM()], such as
#'     `validate`, `reliability` or `scan_polarity`.
#'
#' @return `data.frame` with one row per analysis and columns `"analysis_id"`,
#'     `"mwtab"`, `"toml"`, `"parquet"` and `"mztab"` with the paths to the
#'     files.
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "example_mwtab.json", package = "mwtab2mztabm")
#' res <- mwtabToMzTabM(f, outdir = tempdir(), validate = FALSE)
#' res
#'
#' \donttest{
#' ## Download and convert study ST003934 from the Metabolomics Workbench
#' res <- mwtabToMzTabM("ST003934", outdir = tempdir())
#' }
mwtabToMzTabM <- function(x, outdir = ".", analysis_id = NULL, ...) {
    json <- if (is.character(x) && length(x) == 1L && file.exists(x)) x
            else downloadMwTab(x, destdir = outdir)
    mw <- readMwTab(json)
    if (is.null(analysis_id))
        analysis_id <- names(mw)
    if (!all(analysis_id %in% names(mw)))
        stop("Analysis ", paste(setdiff(analysis_id, names(mw)),
                                collapse = ", "), " not found. Available ",
             "analyses: ", paste(names(mw), collapse = ", "), call. = FALSE)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    res <- lapply(analysis_id, function(an) {
        a <- mw[[an]]
        wb <- a[["METABOLOMICS WORKBENCH"]]
        id <- .mztab_id(wb)
        if (is.na(.first(wb[["ANALYSIS_ID"]])))
            id <- paste0(id, "_", an)
        out <- paste0(file.path(outdir, id), c(".toml", ".parquet", ".mztab"))
        message("Converting analysis ", an)
        withCallingHandlers({
            mwtabToToml(a, out[1L])
            mwtabToParquet(a, out[2L])
            assembleMzTabM(out[1L], out[2L], file = out[3L], ...)
        }, error = function(e) stop("Analysis ", an, ": ",
                                    conditionMessage(e), call. = FALSE))
        data.frame(analysis_id = an, mwtab = json, toml = out[1L],
                   parquet = out[2L], mztab = out[3L])
    })
    do.call(rbind, res)
}

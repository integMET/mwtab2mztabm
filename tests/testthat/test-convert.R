test_that("assembleMzTabM combines the TOML and Parquet files", {
    toml <- example_toml()
    pq <- example_parquet()
    mzt <- assembleMzTabM(toml, pq)
    expect_s4_class(mzt, "MzTabM")
    mtd <- RmzTabM::mtd(mzt)
    expect_equal(mtd_value(mtd, "database[1]-prefix"), "HMDB")
    expect_equal(mtd_value(mtd, "database[2]-prefix"), "KEGG")
    expect_equal(mtd_value(mtd, "mzTab-profile"), "M+S")

    out <- tempfile(fileext = ".mztab")
    assembleMzTabM(toml, pq, file = out, validate = FALSE)
    lines <- readLines(out)
    expect_true(any(startsWith(lines, "MTD\tmzTab-ID\tST999999_AN999998")))
    expect_equal(sum(startsWith(lines, "SMH\t")), 1L)
    expect_equal(sum(startsWith(lines, "SML\t")), 3L)
    ## no tabs within values: all SML rows have as many cells as the header
    ncells <- lengths(strsplit(lines[startsWith(lines, "SM")], "\t"))
    expect_true(all(ncells == ncells[1L]))
})

## mzTab-M file of the example with the value of a cell of the SML section
## replaced, if `value` is provided.
example_mztab <- function(col = NULL, value = NULL) {
    out <- tempfile(fileext = ".mztab")
    assembleMzTabM(example_toml(), example_parquet(), file = out,
                   validate = FALSE)
    if (!is.null(value)) {
        lines <- readLines(out)
        cells <- strsplit(lines[startsWith(lines, "SMH\t")], "\t")[[1L]]
        i <- which(startsWith(lines, "SML\t"))[1L]
        row <- strsplit(lines[i], "\t")[[1L]]
        row[cells == col] <- value
        lines[i] <- paste(row, collapse = "\t")
        writeLines(lines, out)
    }
    out
}

test_that("checkStudyVariableAbundance confirms the mean of the assays", {
    out <- example_mztab()
    res <- checkStudyVariableAbundance(out)
    expect_equal(nrow(res), 0L)
    expect_equal(colnames(res),
                 c("SML_ID", "study_variable", "reported", "expected"))
    expect_true(mwtab2mztabm:::.validate_study_variable_abundance(out))
})

test_that("checkStudyVariableAbundance reports differences", {
    ## study variable 1 (Genotype WT): assays 1 and 3, 1000 and 900.5
    out <- example_mztab("abundance_study_variable[1]", "951")
    res <- checkStudyVariableAbundance(out)
    expect_equal(res, data.frame(SML_ID = 1L,
                                 study_variable = "study_variable[1]",
                                 reported = 951, expected = 950.25))
    expect_error(mwtab2mztabm:::.validate_study_variable_abundance(out),
                 "1 abundance_study_variable value")
    expect_equal(nrow(checkStudyVariableAbundance(out, tolerance = 0.01)), 0L)

    ## an assay value changed instead of the study variable value
    out <- example_mztab("abundance_assay[1]", "2000")
    res <- checkStudyVariableAbundance(out)
    expect_true("study_variable[1]" %in% res$study_variable)
    expect_equal(res$expected[res$study_variable == "study_variable[1]"],
                 1450.25)

    ## a reported value that is missing
    out <- example_mztab("abundance_study_variable[1]", "null")
    res <- checkStudyVariableAbundance(out)
    expect_equal(res$study_variable, "study_variable[1]")
    expect_true(is.na(res$reported))
    expect_equal(res$expected, 950.25)
})

test_that("checkStudyVariableAbundance handles missing assay values", {
    out <- example_mztab()
    ## glucose has no value in assay 2: its means are calculated from the
    ## available values, but expected to be missing with `na.rm = FALSE`.
    expect_equal(nrow(checkStudyVariableAbundance(out)), 0L)
    res <- checkStudyVariableAbundance(out, na.rm = FALSE)
    expect_true(nrow(res) > 0L)
    expect_true(all(res$SML_ID == 2L))
    expect_false(anyNA(res$reported))
    expect_true(all(is.na(res$expected)))
})

test_that("checkStudyVariableAbundance fails for incomplete files", {
    f <- tempfile()
    writeLines("MTD\tmzTab-version\t2.1.0-M", f)
    expect_error(checkStudyVariableAbundance(f), "No SML lines")

    out <- example_mztab()
    lines <- readLines(out)
    lines <- sub("\tabundance_study_variable[1]\t", "\tother\t", lines,
                 fixed = TRUE)
    writeLines(lines, out)
    expect_error(checkStudyVariableAbundance(out), "not found in the SML")
})

test_that("mwtabToMzTabM converts all analyses", {
    d <- tempfile()
    res <- suppressMessages(suppressWarnings(
        mwtabToMzTabM(example_file, outdir = d, validate = FALSE)))
    expect_equal(res$analysis_id, c("AN999998", "AN999999"))
    expect_true(all(file.exists(c(res$toml, res$parquet, res$mztab))))
    expect_equal(basename(res$mztab),
                 c("ST999999_AN999998.mztab", "ST999999_AN999999.mztab"))

    res <- suppressMessages(mwtabToMzTabM(
        example_file, outdir = d, analysis_id = "AN999999", validate = FALSE))
    expect_equal(nrow(res), 1L)
    expect_error(mwtabToMzTabM(example_file, outdir = d, analysis_id = "X"),
                 "not found")
})

test_that("mwtabToMzTabM names files with recovered IDs", {
    mw <- unname(readMwTab(example_file))
    modes <- vapply(mw, function(a) a[["MS"]][["ION_MODE"]], "")
    mw <- lapply(mw, function(a) {
        a[["METABOLOMICS WORKBENCH"]][c("STUDY_ID", "ANALYSIS_ID")] <- NULL
        a
    })
    d <- tempfile()
    dir.create(d)
    f <- file.path(d, "ST000465_mwtab.json")
    jsonlite::write_json(mw, f, auto_unbox = TRUE)
    local_mocked_bindings(.study_analyses = function(id) list(
        list(analysis_id = "AN000726", ion_mode = modes[1L]),
        list(analysis_id = "AN000727", ion_mode = modes[2L])))
    res <- suppressMessages(suppressWarnings(
        mwtabToMzTabM(f, outdir = d, validate = FALSE)))
    expect_equal(res$analysis_id, c("AN000726", "AN000727"))
    expect_equal(basename(res$mztab), c("ST000465_AN000726.mztab",
                                        "ST000465_AN000727.mztab"))
    expect_true(all(file.exists(c(res$toml, res$parquet, res$mztab))))
    lines <- readLines(res$mztab[2L])
    expect_true("MTD\tmzTab-ID\tST000465_AN000727" %in% lines)
})

test_that("mwtabToMzTabM does not overwrite analyses without IDs", {
    mw <- unname(readMwTab(example_file))
    mw <- lapply(mw, function(a) {
        a[["METABOLOMICS WORKBENCH"]][c("STUDY_ID", "ANALYSIS_ID")] <- NULL
        a
    })
    d <- tempfile()
    dir.create(d)
    f <- file.path(d, "local.json")
    jsonlite::write_json(mw, f, auto_unbox = TRUE)
    res <- suppressMessages(suppressWarnings(
        mwtabToMzTabM(f, outdir = d, validate = FALSE)))
    expect_equal(basename(res$mztab), c("mwtab_analysis_1.mztab",
                                        "mwtab_analysis_2.mztab"))
    expect_true(all(file.exists(res$mztab)))
})

test_that("assembleMzTabM keeps non-ASCII text in the C locale", {
    skip_on_os("windows")
    title <- paste0("10 ", intToUtf8(0xb5), "M")
    key <- paste0("BMI [kg/m", intToUtf8(0xb2), "]")
    a <- readMwTab(example_file)[["AN999998"]]
    a$STUDY$STUDY_TITLE <- title
    for (i in seq_along(a$SUBJECT_SAMPLE_FACTORS))
        a$SUBJECT_SAMPLE_FACTORS[[i]]$`Additional sample data`[[key]] <- "27"
    toml <- tempfile(fileext = ".toml")
    pq <- tempfile(fileext = ".parquet")
    out <- tempfile(fileext = ".mztab")
    mwtabToToml(a, toml)
    suppressWarnings(mwtabToParquet(a, pq))

    old <- Sys.getlocale("LC_CTYPE")
    on.exit(Sys.setlocale("LC_CTYPE", old), add = TRUE)
    Sys.setlocale("LC_CTYPE", "C")
    expect_no_warning(assembleMzTabM(toml, pq, file = out, validate = FALSE))
    expect_equal(Sys.getlocale("LC_CTYPE"), "C")

    bytes <- readBin(out, "raw", file.size(out))
    has <- function(x)
        length(grepRaw(charToRaw(enc2utf8(x)), bytes, fixed = TRUE)) > 0L
    expect_true(has(paste0("title\t", title)))
    expect_true(has(paste0("[,,\"", key, "\", 27]")))
    expect_false(has("<U+"))
})

test_that(".with_utf8_ctype restores the locale and warns without UTF-8", {
    skip_on_os("windows")
    old <- Sys.getlocale("LC_CTYPE")
    on.exit(Sys.setlocale("LC_CTYPE", old), add = TRUE)
    Sys.setlocale("LC_CTYPE", "C")

    utf8 <- mwtab2mztabm:::.with_utf8_ctype(l10n_info()$`UTF-8`)
    expect_true(utf8)
    expect_equal(Sys.getlocale("LC_CTYPE"), "C")

    ## also after an error
    expect_error(mwtab2mztabm:::.with_utf8_ctype(stop("failed")), "failed")
    expect_equal(Sys.getlocale("LC_CTYPE"), "C")

    local_mocked_bindings(.UTF8_LOCALES = "not-a-locale")
    expect_warning(res <- mwtab2mztabm:::.with_utf8_ctype(1L),
                   "No UTF-8 locale")
    expect_equal(res, 1L)
    expect_equal(Sys.getlocale("LC_CTYPE"), "C")
})

test_that("converted files pass the jmzTab-M validator", {
    skip_on_cran()
    skip_if_offline("github.com")
    d <- tempfile()
    expect_message(
        res <- suppressWarnings(mwtabToMzTabM(example_file, outdir = d)),
        "Validation finished successfully")
    expect_true(all(file.exists(res$mztab)))
})

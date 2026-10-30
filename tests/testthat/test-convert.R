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

test_that("converted files pass the jmzTab-M validator", {
    skip_on_cran()
    skip_if_offline("github.com")
    d <- tempfile()
    expect_message(
        res <- suppressWarnings(mwtabToMzTabM(example_file, outdir = d)),
        "Validation finished successfully")
    expect_true(all(file.exists(res$mztab)))
})

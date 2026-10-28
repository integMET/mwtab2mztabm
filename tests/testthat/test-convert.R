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

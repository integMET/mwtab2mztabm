test_that("mwtabToParquet exports the quantification table", {
    f <- tempfile(fileext = ".parquet")
    expect_warning(mwtabToParquet(example_file, f, analysis_id = "AN999998"),
                   "1 non-numeric")
    x <- readMwTabParquet(f)
    expect_equal(attr(x, "sample_columns"), c("S1", "S2", "S3", "QC1", "S4"))
    expect_equal(attr(x, "annotation_columns"),
                 c("HMDB", "KEGG ID", "Adduct", "quantified m/z"))
    expect_equal(colnames(x), c("Metabolite", attr(x, "annotation_columns"),
                                attr(x, "sample_columns")))
    expect_equal(x$Metabolite, c("caffeine", "glucose", "unknown_123"))
    expect_true(all(vapply(x[attr(x, "sample_columns")], is.double, NA)))
    expect_equal(x$S3, c(900.5, 450, 11))
    ## "" and "<LOD" are missing values
    expect_equal(x$S2[2L], NA_real_)
    expect_equal(x$S4[2L], NA_real_)
    expect_equal(x$`KEGG ID`, c("C07481", "C00031", ""))
    meta <- attr(x, "metadata")
    expect_equal(meta[["mwtab.units"]], "Peak area")
    expect_equal(meta[["mwtab.analysis_id"]], "AN999998")
})

test_that("mwtabToParquet supports the metabolite_name key", {
    x <- readMwTabParquet(example_parquet("AN999999"))
    expect_equal(attr(x, "annotation_columns"), "pubchem_id")
    expect_equal(x$pubchem_id, c("5950", "612"))
})

test_that("mwtabToParquet fails without metabolite data", {
    a <- readMwTab(example_file)[[1L]]
    a$MS_METABOLITE_DATA <- NULL
    expect_error(mwtabToParquet(a, tempfile()), "does not contain")
})

test_that("annotations are matched by name", {
    meta <- list(list(Metabolite = "b", id = "2"),
                 list(Metabolite = "a", id = "1"))
    res <- mwtab2mztabm:::.metabolite_annotations(meta, c("a", "b"))
    expect_equal(res$id, c("1", "2"))
    expect_warning(
        res <- mwtab2mztabm:::.metabolite_annotations(meta, c("a", "c")),
        "1 metabolite")
    expect_equal(res$id, c("1", NA))

    ## names written differently in the quantification table
    meta <- list(list(Metabolite = "2'-DEOXYCYTIDINE", id = "1"),
                 list(Metabolite = "CITRATE", id = "2"),
                 list(Metabolite = "CITRATE_2", id = "3"),
                 list(Metabolite = "1,3-DIONE", id = "4"))
    expect_warning(expect_message(
        res <- mwtab2mztabm:::.metabolite_annotations(
            meta, c("2^-DEOXYCYTIDINE", "CITRATE", "CITRATE_1", "1_3-DIONE",
                    "1,3-DIONE")),
        "2 metabolite"), "1 metabolite")
    ## "1_3-DIONE" has no unused annotation left
    expect_equal(res$id, c("1", "2", "3", NA, "4"))
})

test_that("readMwTabParquet works without metadata", {
    f <- tempfile(fileext = ".parquet")
    nanoparquet::write_parquet(
        data.frame(Metabolite = c("a", "b"), x = c(1, 2), y = c(3, 4)), f)
    x <- readMwTabParquet(f)
    expect_equal(attr(x, "sample_columns"), c("x", "y"))
    expect_error(readMwTabParquet(f, sample_columns = "z"), "not found")
})

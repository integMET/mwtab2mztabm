test_that("mwtabToToml exports the metadata", {
    f <- example_toml()
    x <- readMwTabToml(f)
    expect_true(all(c("METABOLOMICS WORKBENCH", "PROJECT", "STUDY", "SUBJECT",
                      "COLLECTION", "MS", "SUBJECT_SAMPLE_FACTORS",
                      "MS_METABOLITE_DATA") %in% names(x)))
    ## Only the units of the metabolite data are exported
    expect_equal(x$MS_METABOLITE_DATA, list(Units = "Peak area"))
    expect_equal(x$STUDY$STUDY_SUMMARY,
                 "A synthetic study with two genotypes,\ttwo doses\nand a QC sample.")

    ssf <- x$SUBJECT_SAMPLE_FACTORS
    expect_length(ssf, 5L)
    expect_equal(vapply(ssf, `[[`, "", "Sample ID"),
                 c("S1", "S2", "S3", "QC1", "S4"))
    ## The original order of the factors is restored
    expect_equal(names(ssf[[1L]]$Factors),
                 c("Genotype", "Dose [mg]", "Treatment"))
    expect_equal(ssf[[1L]]$Factors$Treatment, "Control, untreated")
    expect_equal(names(ssf[[1L]][["Additional sample data"]]),
                 c("RAW_FILE_NAME", "injection_order"))
})

test_that("mwtabToToml handles analyses without additional sample data", {
    x <- readMwTabToml(example_toml("AN999999"))
    expect_equal(x$SUBJECT_SAMPLE_FACTORS[[2L]]$Factors, list(Genotype = "KO"))
    expect_length(x$SUBJECT_SAMPLE_FACTORS[[2L]][["Additional sample data"]],
                  0L)
})

test_that("TOML serialization escapes strings and keys", {
    vals <- list(
        a = "quote \" and backslash \\",
        `key with space` = "tab\tnewline\ncarriage\rbell\a",
        `[bracket]` = "unicode ℃",
        arr = c("x", "y"),
        empty = list(),
        num = 1.5, int = 2L, lgl = TRUE,
        tbl = list(b = "1", sub = list(c = "2")),
        aot = list(list(d = "3"), list(d = "4", e = list(f = "5"))))
    f <- tempfile(fileext = ".toml")
    writeLines(enc2utf8(mwtab2mztabm:::.toml_table(list(root = vals))), f,
               useBytes = TRUE)
    res <- unclass(RcppTOML::parseTOML(f, escape = FALSE))$root
    expect_equal(res$a, vals$a)
    expect_equal(res[["key with space"]], vals[["key with space"]])
    expect_equal(res[["[bracket]"]], vals[["[bracket]"]])
    expect_equal(res$arr, c("x", "y"))
    expect_equal(res$num, 1.5)
    expect_equal(res$int, 2L)
    expect_true(res$lgl)
    expect_equal(res$tbl$sub$c, "2")
    expect_equal(res$aot[[2L]]$e$f, "5")
})

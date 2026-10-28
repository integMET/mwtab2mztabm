test_that("mwtabUrl resolves IDs and study page URLs", {
    rest <- "https://www.metabolomicsworkbench.org/rest/study"
    expect_equal(mwtabUrl("ST003934"),
                 paste0(rest, "/study_id/ST003934/mwtab/json"))
    expect_equal(mwtabUrl("an006459"),
                 paste0(rest, "/analysis_id/AN006459/mwtab/json"))
    expect_equal(
        mwtabUrl(paste0("https://www.metabolomicsworkbench.org/data/",
                        "DRCCMetadata.php?Mode=Study&StudyID=ST003934",
                        "&StudyType=MS&ResultType=1")),
        paste0(rest, "/study_id/ST003934/mwtab/json"))
    expect_equal(
        mwtabUrl(paste0("https://www.metabolomicsworkbench.org/data/",
                        "study_textformat_view.php?JSON=YES&",
                        "STUDY_ID=ST003934&ANALYSIS_ID=AN006459")),
        paste0(rest, "/analysis_id/AN006459/mwtab/json"))
    u <- paste0(rest, "/study_id/ST000001/mwtab/json")
    expect_equal(mwtabUrl(u), u)
    expect_equal(mwtabUrl("https://example.org/a.json"),
                 "https://example.org/a.json")
    expect_error(mwtabUrl("foo"), "neither")
    expect_error(mwtabUrl(c("ST000001", "ST000002")), "single")
})

test_that(".is_json_file detects HTML responses", {
    f <- tempfile()
    writeLines("<!DOCTYPE html><html></html>", f)
    expect_false(mwtab2mztabm:::.is_json_file(f))
    writeLines("  [\n{\"a\": 1}]", f)
    expect_true(mwtab2mztabm:::.is_json_file(f))
})

test_that("readMwTab reads all analyses", {
    mw <- readMwTab(example_file)
    expect_equal(names(mw), c("AN999998", "AN999999"))
    expect_true(all(c("STUDY", "SUBJECT_SAMPLE_FACTORS",
                      "MS_METABOLITE_DATA") %in% names(mw[[1L]])))

    ## A single analysis (not within an array) is supported too
    f <- tempfile(fileext = ".json")
    jsonlite::write_json(mw[[2L]], f, auto_unbox = TRUE)
    expect_equal(names(readMwTab(f)), "AN999999")

    writeLines("{\"a\": 1}", f)
    expect_error(readMwTab(f), "not an mwTab")
})

test_that(".mwtab_analysis selects analyses", {
    expect_error(mwtab2mztabm:::.mwtab_analysis(example_file), "2 analyses")
    a <- mwtab2mztabm:::.mwtab_analysis(example_file, "AN999999")
    expect_equal(a[["METABOLOMICS WORKBENCH"]]$ANALYSIS_ID, "AN999999")
    expect_identical(mwtab2mztabm:::.mwtab_analysis(a), a)
    expect_error(mwtab2mztabm:::.mwtab_analysis(example_file, "AN1"),
                 "not found")
})

test_that(".key_values parses factor strings", {
    res <- mwtab2mztabm:::.key_values("Genotype:WT | Treatment:Control")
    expect_equal(res, list(Genotype = "WT", Treatment = "Control"))
    expect_equal(mwtab2mztabm:::.key_values(NULL), setNames(list(),
                                                            character()))
})

test_that("downloadMwTab downloads mwTab JSON", {
    skip_on_cran()
    skip_if_offline("www.metabolomicsworkbench.org")
    d <- tempfile()
    f <- downloadMwTab("AN006459", destdir = d)
    expect_equal(basename(f), "AN006459_mwtab.json")
    expect_equal(names(readMwTab(f)), "AN006459")
})

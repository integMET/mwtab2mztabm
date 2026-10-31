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

## mwTab JSON file of the example without the study and analysis IDs.
example_without_ids <- function(file) {
    mw <- unname(readMwTab(example_file))
    mw <- lapply(mw, function(a) {
        a[["METABOLOMICS WORKBENCH"]][c("STUDY_ID", "ANALYSIS_ID")] <- NULL
        a
    })
    jsonlite::write_json(mw, file, auto_unbox = TRUE)
    file
}

## Analyses of the example as returned by the REST service.
example_rest <- function(ion_modes = c("POSITIVE", "NEGATIVE")) {
    lapply(seq_along(ion_modes), function(i)
        list(analysis_id = paste0("AN00000", i), ion_mode = ion_modes[i]))
}

test_that("readMwTab recovers IDs for files named with a study ID", {
    modes <- vapply(readMwTab(example_file), function(a)
        a[["MS"]][["ION_MODE"]], "", USE.NAMES = FALSE)
    f <- example_without_ids(file.path(tempdir(), "ST000465_mwtab.json"))
    local_mocked_bindings(.study_analyses = function(id) {
        expect_equal(id, "ST000465")
        example_rest(modes)
    })
    mw <- readMwTab(f)
    expect_equal(names(mw), c("AN000001", "AN000002"))
    expect_equal(mw[[2L]][["METABOLOMICS WORKBENCH"]]$STUDY_ID, "ST000465")
    expect_equal(mw[[2L]][["METABOLOMICS WORKBENCH"]]$ANALYSIS_ID, "AN000002")
    expect_equal(names(mw[[1L]][["METABOLOMICS WORKBENCH"]])[1:2],
                 c("STUDY_ID", "ANALYSIS_ID"))
})

test_that("readMwTab keeps the study ID if analyses can not be matched", {
    f <- example_without_ids(file.path(tempdir(), "ST000465_mwtab.json"))
    modes <- list(
        count = "POSITIVE",
        ion_mode = c("NOT", "MATCHING"))
    for (m in modes) {
        local_mocked_bindings(.study_analyses = function(id) example_rest(m))
        expect_warning(mw <- readMwTab(f), "could not be determined")
        expect_equal(names(mw), c("analysis_1", "analysis_2"))
        expect_equal(mw[[1L]][["METABOLOMICS WORKBENCH"]]$STUDY_ID, "ST000465")
        expect_null(mw[[1L]][["METABOLOMICS WORKBENCH"]]$ANALYSIS_ID)
    }
    local_mocked_bindings(.study_analyses = function(id) stop("offline"))
    expect_warning(mw <- readMwTab(f), "could not be determined")
    expect_equal(names(mw), c("analysis_1", "analysis_2"))
})

test_that("readMwTab uses the analysis ID of the file name", {
    mw <- readMwTab(example_file)
    f <- file.path(tempdir(), "AN000726_mwtab.json")
    a <- mw[[1L]]
    a[["METABOLOMICS WORKBENCH"]][c("STUDY_ID", "ANALYSIS_ID")] <- NULL
    jsonlite::write_json(a, f, auto_unbox = TRUE)
    res <- readMwTab(f)
    expect_equal(names(res), "AN000726")
    expect_null(res[[1L]][["METABOLOMICS WORKBENCH"]]$STUDY_ID)
})

test_that("readMwTab warns if IDs can not be determined from the file name", {
    f <- example_without_ids(tempfile(fileext = ".json"))
    expect_warning(mw <- readMwTab(f), "does not contain the study and analysis")
    expect_equal(names(mw), c("analysis_1", "analysis_2"))
})

test_that("readMwTab does not look up IDs of complete files", {
    local_mocked_bindings(.study_analyses = function(id) stop("not expected"))
    expect_equal(names(readMwTab(example_file)), c("AN999998", "AN999999"))
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

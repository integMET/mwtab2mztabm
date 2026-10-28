test_that("mtdFromToml creates the MTD section", {
    mtd <- mtdFromToml(example_toml())
    expect_equal(mtd_value(mtd, "mzTab-ID"), "ST999999_AN999998")
    expect_equal(mtd_value(mtd, "mzTab-profile"), "M+S")
    expect_equal(mtd_value(mtd, "title"), "Synthetic example study: LC-MS")
    ## tabs and newlines are removed
    expect_equal(mtd_value(mtd, "description"),
                 "A synthetic study with two genotypes, two doses and a QC sample.")
    expect_equal(mtd_value(mtd, "contact[1]-name"), "Jane Doe")
    expect_equal(mtd_value(mtd, "contact[1]-email"), "jane.doe@example.org")
    expect_match(mtd_value(mtd, "uri[1]"), "StudyID=ST999999$")
    expect_equal(mtd_value(mtd, "instrument[1]-source"),
                 "[MS, MS:1000073, electrospray ionization, ]")
    expect_equal(mtd_value(mtd, "instrument[1]-analyzer[1]"),
                 "[MS, MS:1000484, orbitrap, ]")
    expect_equal(mtd_value(mtd, "small_molecule-quantification_unit"),
                 "[,, Peak area, ]")
    expect_equal(mtd_value(mtd, "quantification_method"),
                 "[MS, MS:1001834, LC-MS label-free quantitation analysis, ]")
    expect_equal(mtd_value(mtd, "custom[5]"),
                 "[,, chromatography column, \"C18 (100mm x 2.1mm, 1.7um)\"]")
    expect_false(any(mtd[, 1L] == "small_molecule_feature-quantification_unit"))

    ## samples, MS runs and assays
    expect_equal(mtd_value(mtd, "sample[4]"), "QC1")
    expect_equal(mtd_value(mtd, "sample[1]-species[1]"),
                 "[NCBITaxon, NCBITaxon:10090, Mus musculus, ]")
    expect_equal(mtd_value(mtd, "sample[1]-custom[2]"),
                 "[,,RAW_FILE_NAME, s1.mzML]")
    expect_equal(mtd_value(mtd, "ms_run[2]-location"), "s2.mzML")
    expect_equal(mtd_value(mtd, "ms_run[2]-scan_polarity[1]"),
                 "[MS, MS:1000130, positive scan, ]")
    expect_equal(mtd_value(mtd, "assay[5]"), "S4")
    expect_equal(mtd_value(mtd, "assay[5]-ms_run_ref"), "ms_run[5]")
    expect_true(any(mtd_value(mtd, "cv[4]-label") == "NCBITaxon"))

    ## study variables
    expect_equal(mtd_value(mtd, "study_variable_group[2]"),
                 "[,,\"Dose [mg]\",]")
    expect_equal(mtd_value(mtd, "study_variable_group[2]-type"),
                 "[STATO, STATO:0000251, continuous variable, ]")
    expect_equal(mtd_value(mtd, "study_variable_group[1]-type"),
                 "[STATO, STATO:0000252, categorical variable, ]")
    expect_equal(mtd_value(mtd, "study_variable[1]"), "WT")
    expect_equal(mtd_value(mtd, "study_variable[1]-assay_refs"),
                 "assay[1]|assay[3]")
    ## QC sample without dose
    expect_equal(mtd_value(mtd, "study_variable[6]"), "NA")
    expect_equal(mtd_value(mtd, "study_variable[6]-assay_refs"), "assay[4]")
    expect_false(anyNA(mtd))
})

test_that("mtdFromToml orders and subsets assays", {
    expect_message(
        mtd <- mtdFromToml(example_toml(), assays = c("S4", "S1")),
        "3 sample")
    expect_equal(mtd_value(mtd, "assay[1]"), "S4")
    expect_equal(mtd_value(mtd, "ms_run[1]-location"), "s4.mzML")
    expect_false(any(grepl("^assay\\[3\\]", mtd[, 1L])))

    expect_warning(expect_message(
        mtd <- mtdFromToml(example_toml(), assays = c("S1", "X")),
        "4 sample"), "No sample metadata")
    expect_equal(mtd_value(mtd, "assay[2]"), "X")
    expect_equal(mtd_value(mtd, "ms_run[2]-location"), "X")
})

test_that("mtdFromToml uses provided databases and scan polarity", {
    dbs <- data.frame(name = "[,, HMDB, ]", prefix = "HMDB",
                      version = "Unknown", uri = "https://hmdb.ca")
    mtd <- mtdFromToml(example_toml(), databases = dbs,
                       scan_polarity = "negative")
    expect_equal(mtd_value(mtd, "database[1]-prefix"), "HMDB")
    expect_equal(mtd_value(mtd, "ms_run[1]-scan_polarity[1]"),
                 "[MS, MS:1000129, negative scan, ]")
    expect_error(mtdFromToml(example_toml(), databases = dbs[, 1:2]),
                 "needs columns")
})

test_that("mtdFromToml maps GC-MS analyses", {
    mtd <- mtdFromToml(example_toml("AN999999"))
    expect_equal(mtd_value(mtd, "instrument[1]-source"),
                 "[MS, MS:1000389, electron ionization, ]")
    expect_equal(mtd_value(mtd, "instrument[1]-analyzer[1]"),
                 "[MS, MS:1000081, quadrupole, ]")
    expect_equal(mtd_value(mtd, "quantification_method"),
                 "[MS, MS:1001833, quantitation analysis summary, ]")
    expect_equal(mtd_value(mtd, "small_molecule-quantification_unit"),
                 "[PRIDE, PRIDE:0000330, Arbitrary quantification unit, ]")
    ## no raw file name available: sample ID is used as location
    expect_equal(mtd_value(mtd, "ms_run[1]-location"), "S1")
})

test_that("mtdFromToml maps hybrid analyzers and skips placeholders", {
    x <- readMwTabToml(example_toml())
    x$MS$INSTRUMENT_TYPE <- "QTOF"
    x$CHROMATOGRAPHY$CHROMATOGRAPHY_TYPE <- "-"
    x$CHROMATOGRAPHY$COLUMN_NAME <- "-"
    mtd <- mtdFromToml(x)
    expect_equal(mtd_value(mtd, "instrument[1]-analyzer[1]"),
                 "[MS, MS:1000081, quadrupole, ]")
    expect_equal(mtd_value(mtd, "instrument[1]-analyzer[2]"),
                 "[MS, MS:1000084, time-of-flight, ]")
    expect_false(any(grepl("chromatography (type|column)", mtd[, 2L])))
    expect_equal(mtd_value(mtd, "quantification_method"),
                 "[MS, MS:1001833, quantitation analysis summary, ]")

    an <- function(x) mwtab2mztabm:::.map_cv(
        x, mwtab2mztabm:::.ANALYZER, "MS:1000443", "mass analyzer type",
        all = TRUE)
    expect_equal(an("Q Exactive"), c("[MS, MS:1000081, quadrupole, ]",
                                     "[MS, MS:1000484, orbitrap, ]"))
    expect_equal(an("Triple quadrupole"), "[MS, MS:1000081, quadrupole, ]")
    expect_equal(an("GC-TOF"), "[MS, MS:1000084, time-of-flight, ]")
    expect_equal(an("-"), character())
})

test_that("mtdFromToml requires a known scan polarity", {
    x <- readMwTabToml(example_toml())
    x$MS$ION_MODE <- "UNSPECIFIED"
    expect_error(mtdFromToml(x), "scan polarity")
    mtd <- mtdFromToml(x, scan_polarity = "positive")
    expect_equal(mtd_value(mtd, "ms_run[1]-scan_polarity[1]"),
                 "[MS, MS:1000130, positive scan, ]")
})

test_that("CV parameter values are quoted if needed", {
    expect_equal(mwtab2mztabm:::.cv_quote(c("a", "a, b", "[x]", "a\"b")),
                 c("a", "\"a, b\"", "\"[x]\"", "a'b"))
    expect_equal(mwtab2mztabm:::.user_param("unit", "mg, dry"),
                 "[,, unit, \"mg, dry\"]")
})

test_that("mtdFromToml works without species and sample type", {
    x <- readMwTabToml(example_toml())
    x$SUBJECT <- NULL
    x$COLLECTION <- NULL
    mtd <- mtdFromToml(x)
    expect_false(any(grepl("species|tissue", mtd[, 1L])))
    expect_false("NCBITaxon" %in% mtd[, 2L])
})

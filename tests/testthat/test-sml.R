test_that("smlFromParquet creates the SML section", {
    pq <- example_parquet()
    sml <- smlFromParquet(pq)
    expect_equal(nrow(sml), 3L)
    expect_equal(sml$database_identifier,
                 c("HMDB:HMDB0001847", "KEGG:C00031", "null"))
    expect_equal(sml$uri, c("https://hmdb.ca/metabolites/HMDB0001847",
                            "https://www.kegg.jp/entry/C00031", "null"))
    ## names only for molecules with database identifiers
    expect_equal(sml$chemical_name, c("caffeine", "glucose", "null"))
    expect_equal(sml$opt_global_metabolite_name,
                 c("caffeine", "glucose", "unknown_123"))
    expect_equal(sml$adduct_ions, c("[M+H]1+", "[M+Na]1+", "null"))
    expect_equal(sml$reliability, rep("2", 3L))
    expect_true(all(c("opt_global_hmdb", "opt_global_kegg_id",
                      "opt_global_adduct", "opt_global_quantified_m_z")
                    %in% colnames(sml)))
    expect_equal(sml[["abundance_assay[3]"]], c(900.5, 450, 11))
    expect_true(is.na(sml[["abundance_assay[2]"]][2L]))
    expect_false(any(grepl("study_variable", colnames(sml))))

    sml <- smlFromParquet(pq, reliability = "3")
    expect_equal(sml$reliability, rep("3", 3L))
})

test_that("smlFromParquet adds study variable columns", {
    tab <- readMwTabParquet(example_parquet())
    mtd <- mtdFromToml(example_toml(), assays = attr(tab, "sample_columns"))
    sml <- smlFromParquet(tab, mtd = mtd)
    expect_true(all(paste0("abundance_study_variable[", 1:9, "]") %in%
                    colnames(sml)))
    ## study variable 1 (Genotype WT): assays 1 and 3
    expect_equal(sml[["abundance_study_variable[1]"]],
                 rowMeans(cbind(c(1000, 500, 10), c(900.5, 450, 11))))
    ## means ignore missing abundances: glucose has no value in assay 2
    abundance <- as.matrix(sml[, paste0("abundance_assay[", 1:3, "]")])
    expect_true(is.na(abundance[2L, 2L]))
    for (sv in grep("^study_variable\\[[0-9]+\\]-assay_refs$", mtd[, 1L],
                    value = TRUE)) {
        refs <- strsplit(mtd_value(mtd, sv), "|", fixed = TRUE)[[1L]]
        col <- paste0("abundance_", sub("-assay_refs$", "", sv))
        expect_equal(sml[[col]], unname(rowMeans(
            as.matrix(sml[, paste0("abundance_", refs), drop = FALSE]),
            na.rm = TRUE)))
    }
    ## mean is only missing if all abundances are missing
    tab2 <- tab
    tab2[1L, attr(tab, "sample_columns")] <- NA_real_
    sml2 <- smlFromParquet(tab2, mtd = mtd)
    expect_true(is.na(sml2[["abundance_study_variable[1]"]][1L]))
    expect_false(is.nan(sml2[["abundance_study_variable[1]"]][1L]))
    expect_equal(sml2[["abundance_study_variable[1]"]][2:3],
                 sml[["abundance_study_variable[1]"]][2:3])

    mtd2 <- suppressMessages(mtdFromToml(example_toml(),
                                         assays = c("S1", "S2")))
    expect_error(smlFromParquet(tab, mtd = mtd2), "2 assays")
})

test_that("smlDatabases lists the referenced databases", {
    pq <- example_parquet()
    dbs <- smlDatabases(pq)
    expect_equal(dbs$prefix, c("HMDB", "KEGG"))
    expect_equal(colnames(dbs), c("name", "prefix", "version", "uri"))
    expect_equal(smlDatabases(example_parquet("AN999999"))$prefix, "PubChem")
})

test_that(".adduct normalizes adducts", {
    expect_equal(mwtab2mztabm:::.adduct(
                     c("[M+H]+", "[M-H]-", "[M+2H]2+", "[2M+Na]+", "M+H", "")),
                 c("[M+H]1+", "[M-H]1-", "[M+2H]2+", "[2M+Na]1+", NA, NA))
})

test_that(".opt_name creates valid column names", {
    expect_equal(mwtab2mztabm:::.opt_name(
                     c("quantified m/z", "KEGG ID", "_x_", "%")),
                 c("quantified_m_z", "kegg_id", "x", "annotation"))
})

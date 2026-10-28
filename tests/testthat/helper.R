example_file <- system.file("extdata", "example_mwtab.json",
                            package = "mwtab2mztabm")

## TOML and Parquet files of an analysis of the example mwTab file.
example_toml <- function(analysis_id = "AN999998") {
    f <- tempfile(fileext = ".toml")
    mwtabToToml(example_file, f, analysis_id = analysis_id)
    f
}

example_parquet <- function(analysis_id = "AN999998") {
    f <- tempfile(fileext = ".parquet")
    suppressWarnings(mwtabToParquet(example_file, f,
                                    analysis_id = analysis_id))
    f
}

mtd_value <- function(mtd, field) {
    unname(mtd[mtd[, 1L] == field, 2L])
}

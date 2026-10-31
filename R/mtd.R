#' Create the mzTab-M metadata (MTD) section from the TOML file
#'
#' `mtdFromToml()` creates the mzTab-M MTD section from the mwTab metadata
#' exported with [mwtabToToml()] using the functions of the *RmzTabM* package.
#' The mwTab content is mapped as follows:
#'
#' - `mzTab-ID`: `"<STUDY_ID>_<ANALYSIS_ID>"`.
#' - `title`, `description`: `STUDY_TITLE` and `STUDY_SUMMARY` (falling back
#'   to the respective `PROJECT` fields).
#' - `contact[1]`: name, institute and email of the `STUDY` (or `PROJECT`).
#' - `publication[1]`: the `DOI` of the `PROJECT`, if available.
#' - `uri[1]`: the Metabolomics Workbench page of the study.
#' - `instrument[1]`: the `INSTRUMENT_NAME`, ionization (`MS_TYPE`) and mass
#'   analyzer (`INSTRUMENT_TYPE`) from section `MS`.
#' - `custom[1-n]`: project and analysis ID and the chromatography type,
#'   instrument and column.
#' - `sample[1-n]`: one sample for each element of `SUBJECT_SAMPLE_FACTORS`,
#'   with the species from section `SUBJECT`, the sample type from section
#'   `COLLECTION` as tissue and the subject ID and additional sample data as
#'   custom parameters.
#' - `ms_run[1-n]`: one MS run per sample, with the location taken from the
#'   raw file name in the additional sample data (key `RAW_FILE_NAME`) and the
#'   scan polarity from the `ION_MODE`. If no raw file name is available, the
#'   sample ID is used as location since the jmzTab-M validator does not
#'   accept `"null"`. The location is percent-encoded to be a valid URI.
#' - `assay[1-n]`: one assay per sample, linking to its sample and MS run.
#' - `study_variable_group[1-n]`, `study_variable[1-n]`: one study variable
#'   group for each sample factor, with one study variable for each of its
#'   distinct values. Factors with only numeric values are reported as
#'   continuous variables, all others as categorical variables.
#' - `small_molecule-quantification_unit`: the `Units` of the metabolite data.
#'
#' @param x `character(1)` with the path to the TOML file created by
#'     [mwtabToToml()] or the result of [readMwTabToml()].
#'
#' @param assays optional `character` with the sample IDs to report as
#'     assays, in the order of the abundance columns of the SML section. By
#'     default, all samples are reported in the order of the TOML file.
#'
#' @param databases optional `data.frame` with columns `"name"`, `"prefix"`,
#'     `"version"` and `"uri"` defining the databases referenced in the
#'     `database_identifier` column of the SML section. If not provided,
#'     `"no database"` is reported.
#'
#' @param scan_polarity optional `character(1)`, either `"positive"` or
#'     `"negative"`. Defaults to the `ION_MODE` of section `MS`.
#'
#' @return Two-column `character` `matrix` with the MTD section.
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "example_mwtab.json", package = "mwtab2mztabm")
#' toml <- tempfile(fileext = ".toml")
#' mwtabToToml(f, toml, analysis_id = "AN999998")
#'
#' mtd <- mtdFromToml(toml)
#' head(mtd, 20)
mtdFromToml <- function(x, assays = NULL, databases = NULL,
                        scan_polarity = NULL) {
    if (is.character(x))
        x <- readMwTabToml(x)
    mw <- x[["METABOLOMICS WORKBENCH"]]
    study <- x[["STUDY"]]
    project <- x[["PROJECT"]]
    ms <- x[["MS"]]
    species <- .species_param(x[["SUBJECT"]])
    cvs <- .cv_table(ncbitaxon = any(grepl("NCBITaxon:", species)))
    dbs <- .database_table(databases)
    mtd <- RmzTabM::mtdSkeleton(
        id = .mztab_id(mw),
        software = .cv_param("MS", "MS:1000799",
                             "custom unreleased software tool",
                             paste("mwtab2mztabm", .pkg_version())),
        quantification_method = .quantification_method(x[["CHROMATOGRAPHY"]]),
        cv_label = cvs$label, cv_full_name = cvs$full_name,
        cv_version = cvs$version, cv_uri = cvs$uri,
        database = dbs$name, database_prefix = dbs$prefix,
        database_version = dbs$version, database_uri = dbs$uri,
        small_molecule_quantification_unit = .quantification_unit(x),
        mztab_profile = "M+S")
    ## Only needed for files with a SMF section.
    mtd <- mtd[mtd[, 1L] != "small_molecule_feature-quantification_unit", ,
               drop = FALSE]
    title <- .first(study$STUDY_TITLE, .first(project$PROJECT_TITLE, ""))
    descr <- .first(study$STUDY_SUMMARY,
                    .first(project$PROJECT_SUMMARY, title))
    mtd <- rbind(mtd,
                 c("title", .clean_text(title)),
                 c("description", .clean_text(descr)),
                 .contact_fields(study, project),
                 .publication_fields(project),
                 .uri_fields(mw),
                 .instrument_fields(ms),
                 .custom_fields(mw, x[["CHROMATOGRAPHY"]]))
    sdata <- .sample_data(x, assays)
    has_instrument <- any(grepl("^instrument\\[1\\]", mtd[, 1L]))
    sdata$..location <- .location_uri(.clean_text(sdata$..location))
    sdata$..scan_polarity <- scan_polarity %||% .scan_polarity(ms)
    if (length(species))
        sdata$..species <- species
    tissue <- .tissue_param(x[["COLLECTION"]])
    if (length(tissue))
        sdata$..tissue <- tissue
    if (has_instrument)
        sdata$..instrument_ref <- 1L
    custom <- grep("^\\.\\.custom_", colnames(sdata), value = TRUE)
    factors <- attr(sdata, "factors")
    ## Columns are mapped explicitly (not with RmzTabM::sampleCols() etc) to
    ## avoid sample factors being picked up by the default column names.
    scols <- c(sample = "..sample", species = "..species",
               tissue = "..tissue")
    smtd <- RmzTabM::mtdFromSampleData(
        sdata,
        sampleCols. = c(scols[scols %in% colnames(sdata)],
                        setNames(custom, attr(sdata, "custom_names"))),
        msRunCols. = c(location = "..location",
                       scan_polarity = "..scan_polarity",
                       instrument_ref = "..instrument_ref"),
        assayCols. = c(assay = "..sample"),
        groups = factors,
        group_description = if (length(factors))
            paste0("mwTab sample factor ", factors) else character())
    ## Study variable for samples with missing factor values.
    smtd[is.na(smtd[, 2L]), 2L] <- "NA"
    RmzTabM::mtdSort(rbind(mtd, smtd))
}

## Sample information as `data.frame`, one row per assay. Columns starting
## with `".."` are fixed; factor columns are named by the (quoted) factor
## names; additional sample data are in columns `"..custom_<i>"`.
.sample_data <- function(x, assays = NULL) {
    ssf <- x[["SUBJECT_SAMPLE_FACTORS"]]
    ids <- vapply(ssf, function(z) .first(z[["Sample ID"]], ""), "")
    if (any(d <- duplicated(ids)))
        warning("Duplicated sample ID(s) ", paste(unique(ids[d]),
                                                    collapse = ", "),
                ". Only the first entry is used.", call. = FALSE)
    if (is.null(assays))
        assays <- ids[!d]
    idx <- match(assays, ids)
    if (anyNA(idx))
        warning("No sample metadata available for ", sum(is.na(idx)),
                " sample(s): ", paste(head(assays[is.na(idx)], 5L),
                                      collapse = ", "), call. = FALSE)
    if (any(dropped <- !ids[!d] %in% assays))
        message(sum(dropped), " sample(s) in the sample metadata have no ",
                "abundance values and are not reported")
    ssf <- ssf[idx]
    ord <- x[["mwtab2mztabm"]]
    fnames <- .names_in_order(ssf, "Factors", ord$factor_order)
    anames <- .names_in_order(ssf, "Additional sample data",
                              ord$additional_data_order)
    getv <- function(element, key) vapply(ssf, function(z) {
        v <- if (is.null(z)) NULL else z[[element]][[key]]
        .first(v)
    }, NA_character_)
    res <- data.frame(..sample = .clean_text(assays))
    raw <- grep("^raw[ _]?file", anames, ignore.case = TRUE, value = TRUE)
    res$..location <- if (length(raw)) getv("Additional sample data", raw[1L])
                      else NA_character_
    ## The jmzTab-M validator does not accept "null" as location.
    nofile <- .is_na_value(res$..location)
    res$..location[nofile] <- res$..sample[nofile]
    custom <- list()
    subj <- vapply(ssf, function(z) .first(z[["Subject ID"]]), NA_character_)
    if (!all(.is_na_value(subj)))
        custom[["Subject ID"]] <- subj
    for (a in anames)
        custom[[a]] <- getv("Additional sample data", a)
    for (i in seq_along(custom)) {
        v <- custom[[i]]
        v[is.na(v)] <- ""
        res[[paste0("..custom_", i)]] <- .cv_quote(v)
    }
    for (f in fnames)
        res[[.cv_quote(f)]] <- .factor_values(getv("Factors", f))
    ## Avoid clashes with the parameter names of RmzTabM::mtdSample().
    cnames <- names(custom)
    clash <- cnames %in% names(formals(RmzTabM::mtdSample))
    cnames[clash] <- paste(cnames[clash], "(mwTab)")
    attr(res, "custom_names") <- .cv_quote(cnames)
    attr(res, "factors") <- .cv_quote(fnames)
    res
}

## The location of an MS run has to be a valid URI, which file names and sample
## IDs often are not (e.g. `"Sample # 1"`). Their reserved characters are thus
## percent-encoded, except for path separators. Of a URI, only characters that
## are not allowed in URIs are encoded, and values that are already
## percent-encoded are not changed.
.location_uri <- function(x) {
    vapply(enc2utf8(x), function(u) {
        if (grepl("%[[:xdigit:]]{2}", u))
            return(u)
        if (grepl("^[[:alpha:]][[:alnum:]+.-]*://", u))
            return(URLencode(u))
        gsub("%2F", "/", URLencode(u, reserved = TRUE), fixed = TRUE)
    }, "", USE.NAMES = FALSE)
}

.names_in_order <- function(ssf, element, ord) {
    nms <- unique(unlist(lapply(ssf, function(z) names(z[[element]])),
                         use.names = FALSE)) %||% character()
    c(intersect(ord, nms), setdiff(nms, ord))
}

.is_na_value <- function(x) {
    is.na(x) | trimws(x) %in% .NA_STRINGS
}

## Factor values: numeric if all non-missing values are numbers, character
## otherwise.
.factor_values <- function(x) {
    x <- .clean_text(x)
    num <- suppressWarnings(as.numeric(x))
    if (all(!is.na(num) | .is_na_value(x)) && any(!is.na(num)))
        return(num)
    x[!nzchar(x)] <- NA_character_
    x
}

.mztab_id <- function(mw) {
    ids <- c(.first(mw$STUDY_ID), .first(mw$ANALYSIS_ID))
    ids <- ids[!is.na(ids)]
    if (length(ids)) paste(ids, collapse = "_") else "mwtab"
}

.scan_polarity <- function(ms) {
    mode <- tolower(.first(ms$ION_MODE, ""))
    if (grepl("^pos", mode)) return("positive")
    if (grepl("^neg", mode)) return("negative")
    stop("Can not determine the scan polarity from ION_MODE \"",
         .first(ms$ION_MODE, ""), "\". Please provide it with parameter ",
         "'scan_polarity'.", call. = FALSE)
}

.species_param <- function(subject) {
    sp <- .first(subject$SUBJECT_SPECIES, "")
    if (!nzchar(sp))
        return(character())
    taxid <- .first(subject$TAXONOMY_ID, "")
    if (grepl("^[0-9]+$", taxid))
        .cv_param("NCBITaxon", paste0("NCBITaxon:", taxid), .cv_quote(sp))
    else .user_param(sp)
}

.tissue_param <- function(collection) {
    st <- .first(collection$SAMPLE_TYPE, "")
    if (nzchar(st)) .user_param(st) else character()
}

## mwTab CHROMATOGRAPHY_TYPE is e.g. "Reversed phase", "HILIC" or "GC".
.quantification_method <- function(chrom) {
    type <- .first(chrom$CHROMATOGRAPHY_TYPE, "")
    if (nzchar(type) && !.is_na_value(type) &&
        !grepl("\\bGC\\b|^gas|\\bCE\\b|electrophoresis|none|direct",
               type, ignore.case = TRUE))
        "[MS, MS:1001834, LC-MS label-free quantitation analysis, ]"
    else "[MS, MS:1001833, quantitation analysis summary, ]"
}

.quantification_unit <- function(x) {
    units <- unlist(lapply(x[.data_section_names(x)], `[[`, "Units"))
    units <- .first(units, "")
    if (.is_na_value(units) ||
        tolower(units) %in% c("not applicable", "none", "na", "n/a"))
        "[PRIDE, PRIDE:0000330, Arbitrary quantification unit, ]"
    else .user_param(units)
}

.contact_fields <- function(study, project) {
    src <- if (!.is_blank(study$LAST_NAME)) study else project
    name <- .clean_text(paste(.first(src$FIRST_NAME, ""),
                              .first(src$LAST_NAME, "")))
    vals <- list(name = name,
                 affiliation = .clean_text(.first(src$INSTITUTE, "")),
                 email = .clean_text(.first(src$EMAIL, "")))
    vals <- vals[nzchar(unlist(vals))]
    if (!length(vals)) return(NULL)
    do.call(RmzTabM::mtdFields, c(vals, field_prefix = "contact"))
}

.publication_fields <- function(project) {
    doi <- sub("^(https?://(dx\\.)?doi\\.org/|doi:)", "",
               .first(project$DOI, ""), ignore.case = TRUE)
    if (grepl("^10\\.", doi)) c("publication[1]", paste0("doi:", doi))
}

.uri_fields <- function(mw) {
    st <- .first(mw$STUDY_ID, "")
    if (nzchar(st))
        c("uri[1]", paste0("https://www.metabolomicsworkbench.org/data/",
                           "DRCCMetadata.php?Mode=Study&StudyID=", st))
}

## Ionization types (mwTab MS_TYPE) and mass analyzers (mwTab
## INSTRUMENT_TYPE) mapped to PSI-MS terms; patterns are matched in order.
.IONIZATION <- rbind(
    c("^ESI|electrospray", "MS:1000073", "electrospray ionization"),
    c("^APCI|atmospheric pressure chemical", "MS:1000070",
      "atmospheric pressure chemical ionization"),
    c("^APPI|photoioni", "MS:1000382", "atmospheric pressure photoionization"),
    c("^MALDI|laser desorption", "MS:1000075",
      "matrix-assisted laser desorption ionization"),
    c("^EI$|electron ioni", "MS:1000389", "electron ionization"),
    c("^CI$|chemical ioni", "MS:1000071", "chemical ionization"))

## All matching analyzers are reported for hybrid instruments (e.g. "QTOF").
.ANALYZER <- rbind(
    c("quadrupole|^q$|qqq|\\bqq?-?tof|qtrap|\\bq[ -]?exactive|triple ?tof",
      "MS:1000081", "quadrupole"),
    c("ion trap|qtrap|^it$|\\bltq", "MS:1000264", "ion trap"),
    c("orbitrap|exactive", "MS:1000484", "orbitrap"),
    c("tof|time.of.flight", "MS:1000084", "time-of-flight"),
    c("ft-?icr|cyclotron", "MS:1000079",
      "fourier transform ion cyclotron resonance"))

.map_cv <- function(x, map, parent_accession, parent_name, all = FALSE) {
    if (.is_na_value(x)) return(character())
    hit <- which(vapply(map[, 1L], grepl, NA, x = x, ignore.case = TRUE))
    if (!all)
        hit <- head(hit, 1L)
    if (length(hit))
        .cv_param("MS", map[hit, 2L], map[hit, 3L])
    else .cv_param("MS", parent_accession, parent_name, x)
}

.instrument_fields <- function(ms) {
    name <- .first(ms$INSTRUMENT_NAME, "")
    if (.is_na_value(name)) return(NULL)
    analyzers <- .map_cv(.first(ms$INSTRUMENT_TYPE, ""), .ANALYZER,
                         "MS:1000443", "mass analyzer type", all = TRUE)
    vals <- c(list(
        name = .cv_param("MS", "MS:1000031", "instrument model", name),
        source = .map_cv(.first(ms$MS_TYPE, ""), .IONIZATION,
                         "MS:1000008", "ionization type")),
        setNames(as.list(analyzers),
                 sprintf("analyzer[%d]", seq_along(analyzers))))
    vals <- vals[lengths(vals) > 0L]
    do.call(RmzTabM::mtdFields, c(vals, field_prefix = "instrument"))
}

.custom_fields <- function(mw, chrom) {
    vals <- c(
        `Metabolomics Workbench project ID` = .first(mw$PROJECT_ID, ""),
        `Metabolomics Workbench analysis ID` = .first(mw$ANALYSIS_ID, ""),
        `chromatography type` = .first(chrom$CHROMATOGRAPHY_TYPE, ""),
        `chromatography instrument` = .first(chrom$INSTRUMENT_NAME, ""),
        `chromatography column` = .first(chrom$COLUMN_NAME, ""))
    vals <- vals[!.is_na_value(vals)]
    if (!length(vals)) return(NULL)
    cbind(paste0("custom[", seq_along(vals), "]"),
          .user_param(names(vals), vals))
}

.cv_table <- function(ncbitaxon = FALSE) {
    f <- formals(RmzTabM::mtdSkeleton)
    res <- data.frame(label = eval(f$cv_label),
                      full_name = eval(f$cv_full_name),
                      version = eval(f$cv_version),
                      uri = eval(f$cv_uri))
    if (ncbitaxon)
        res <- rbind(res, data.frame(
            label = "NCBITaxon", full_name = "NCBI organismal classification",
            version = "unknown",
            uri = "http://purl.obolibrary.org/obo/ncbitaxon.owl"))
    res
}

.database_table <- function(x = NULL) {
    if (is.null(x) || !nrow(x)) {
        f <- formals(RmzTabM::mtdSkeleton)
        return(data.frame(name = eval(f$database),
                          prefix = eval(f$database_prefix),
                          version = eval(f$database_version),
                          uri = eval(f$database_uri)))
    }
    if (!all(c("name", "prefix", "version", "uri") %in% colnames(x)))
        stop("'databases' needs columns \"name\", \"prefix\", \"version\" ",
             "and \"uri\"", call. = FALSE)
    x
}

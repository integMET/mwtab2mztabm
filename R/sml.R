#' Create the mzTab-M small molecule (SML) section from the Parquet file
#'
#' `smlFromParquet()` creates the mzTab-M SML section from the metabolite
#' quantification table exported with [mwtabToParquet()] using
#' [RmzTabM::smlCreate()]. Each metabolite is reported as one small molecule
#' with its abundances in each assay. If the MTD section is provided, the
#' average and variation of the abundances per study variable are added with
#' [RmzTabM::smlAddStudyVariableColumns()].
#'
#' The metabolite annotations are mapped as follows:
#'
#' - `database_identifier`, `uri`: the first available identifier from an
#'   HMDB (`"HMDB:HMDB0000001"`), KEGG (`"KEGG:C00001"`) or PubChem
#'   (`"PubChem:5793"`) annotation column (matched by name, e.g. `"HMDB"`,
#'   `"hmdb_id"`, `"kegg_id"`, `"pubchem_id"`), or `"null"`.
#' - `chemical_name`: the mwTab metabolite name. The mzTab-M validator requires
#'   the number of names to match the number of database identifiers, so
#'   names are only reported for metabolites with a database identifier.
#' - `chemical_formula`, `smiles`, `inchi`: from respective annotation columns,
#'   if available (only for metabolites with a database identifier).
#' - `adduct_ions`: from an `"adduct"` annotation column, if available and in
#'   a valid format (e.g. `"[M+H]1+"`).
#' - `reliability`: the value of parameter `reliability`. mwTab files do not
#'   report the identification confidence of metabolites, but the column must
#'   not be empty.
#' - `opt_global_metabolite_name` and `opt_global_<annotation>`: the mwTab
#'   metabolite name and all annotation columns as optional columns.
#'
#' `smlDatabases()` returns the definition of the databases referenced in the
#' `database_identifier` column, for use in the MTD section (see parameter
#' `databases` of [mtdFromToml()]).
#'
#' @param x `character(1)` with the path to the Parquet file created by
#'     [mwtabToParquet()] or the result of [readMwTabParquet()].
#'
#' @param mtd optional MTD section (see [mtdFromToml()]) used to add the
#'     study variable abundance columns. Its assays have to match the sample
#'     columns of `x`.
#'
#' @param reliability `character(1)` with the identification reliability
#'     reported for all metabolites. Defaults to `"2"` (*putatively annotated
#'     compound*) of the compound identification confidence levels.
#'
#' @return `smlFromParquet()`: `data.frame` with the SML section.
#'     `smlDatabases()`: `data.frame` with columns `"name"`, `"prefix"`,
#'     `"version"` and `"uri"`.
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "example_mwtab.json", package = "mwtab2mztabm")
#' pq <- tempfile(fileext = ".parquet")
#' mwtabToParquet(f, pq, analysis_id = "AN999998")
#'
#' sml <- smlFromParquet(pq)
#' sml[, c("SML_ID", "database_identifier", "chemical_name",
#'         "abundance_assay[1]")]
#' smlDatabases(pq)
smlFromParquet <- function(x, mtd = NULL, reliability = "2") {
    if (is.character(x))
        x <- readMwTabParquet(x)
    samples <- attr(x, "sample_columns")
    ann <- .sml_annotations(x)
    amat <- as.matrix(x[, samples, drop = FALSE])
    storage.mode(amat) <- "double"
    args <- c(.sml_opt_columns(x),
              list(x = amat,
                   database_identifier = ann$id,
                   chemical_name = ann$name,
                   uri = ann$uri,
                   chemical_formula = ann$formula,
                   smiles = ann$smiles,
                   inchi = ann$inchi,
                   adduct_ions = ann$adduct,
                   reliability = rep(as.character(reliability), nrow(x))))
    sml <- do.call(RmzTabM::smlCreate, args)
    if (!is.null(mtd)) {
        n <- sum(grepl("^assay\\[[0-9]+\\]$", mtd[, 1L]))
        if (n != length(samples))
            stop("The MTD section defines ", n, " assays but the ",
                 "quantification table has ", length(samples), " samples",
                 call. = FALSE)
        sml <- RmzTabM::smlAddStudyVariableColumns(sml, mtd)
    }
    sml
}

#' @rdname smlFromParquet
#'
#' @export
smlDatabases <- function(x) {
    if (is.character(x))
        x <- readMwTabParquet(x)
    used <- unique(sub(":.*$", "", .sml_annotations(x)$id))
    .DATABASES[.DATABASES$prefix %in% used,
               c("name", "prefix", "version", "uri")]
}

## Supported databases with the (normalized) names of annotation columns
## providing their identifiers and the pattern of valid identifiers.
.DATABASES <- data.frame(
    prefix = c("HMDB", "KEGG", "PubChem"),
    name = c("[,, HMDB, ]", "[,, KEGG, ]", "[,, PubChem, ]"),
    version = "Unknown",
    uri = c("https://hmdb.ca", "https://www.kegg.jp",
            "https://pubchem.ncbi.nlm.nih.gov"),
    columns = I(list(c("hmdb", "hmdbid", "hmdbids", "hmdbidentifier"),
                     c("kegg", "keggid", "keggids", "keggcompoundid"),
                     c("pubchem", "pubchemid", "pubchemcid", "cid",
                       "pubchemcompoundid"))),
    pattern = c("^HMDB[0-9]+$", "^[CDG][0-9]{5}$", "^[0-9]+$"),
    entry_uri = c("https://hmdb.ca/metabolites/",
                  "https://www.kegg.jp/entry/",
                  "https://pubchem.ncbi.nlm.nih.gov/compound/"))

.norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

## Annotation column (character) whose normalized name is one of `names`.
.annotation_column <- function(x, names) {
    acols <- attr(x, "annotation_columns")
    hit <- acols[.norm_name(acols) %in% names]
    if (length(hit)) .clean_text(x[[hit[1L]]]) else rep("", nrow(x))
}

.sml_annotations <- function(x) {
    n <- nrow(x)
    id <- uri <- rep(NA_character_, n)
    for (i in seq_len(nrow(.DATABASES))) {
        db <- .DATABASES[i, ]
        v <- .annotation_column(x, db$columns[[1L]])
        ## Use the first of possibly multiple IDs and drop a database prefix.
        v <- sub("^[A-Za-z]+:", "", sub("[;,| ].*$", "", v))
        if (db$prefix == "HMDB") v <- toupper(v)
        ok <- is.na(id) & grepl(db$pattern, v)
        id[ok] <- paste0(db$prefix, ":", v[ok])
        uri[ok] <- paste0(db$entry_uri, v[ok])
    }
    has_id <- !is.na(id)
    only_id <- function(v) {
        v[!has_id | !nzchar(v)] <- NA_character_
        v
    }
    list(id = id,
         uri = uri,
         name = only_id(gsub("|", "/", .clean_text(x$Metabolite),
                             fixed = TRUE)),
         formula = only_id(.annotation_column(
             x, c("formula", "chemicalformula", "molecularformula"))),
         smiles = only_id(.annotation_column(
             x, c("smiles", "canonicalsmiles", "isomericsmiles"))),
         inchi = only_id(.annotation_column(x, c("inchi", "standardinchi"))),
         adduct = .adduct(.annotation_column(
             x, c("adduct", "adducts", "adduction", "adductions"))))
}

## Adducts in mzTab-M format, e.g. `"[M+H]1+"`; `NA` if not valid.
.adduct <- function(x) {
    x <- gsub(" ", "", x, fixed = TRUE)
    ok <- grepl("^\\[[^]]*M[^]]*\\][0-9]*[+-]$", x)
    x[!ok] <- NA_character_
    sub("\\]([+-])$", "]1\\1", x)
}

## Optional columns: the metabolite name and all annotation columns.
.sml_opt_columns <- function(x) {
    acols <- attr(x, "annotation_columns")
    res <- lapply(c("Metabolite", acols), function(z) .clean_text(x[[z]]))
    nms <- .opt_name(c("metabolite_name", acols))
    names(res) <- paste0("opt_global_", make.unique(nms, sep = "_"))
    res
}

.opt_name <- function(x) {
    x <- gsub("[^a-z0-9]+", "_", tolower(x))
    x <- gsub("^_+|_+$", "", x)
    x[!nzchar(x)] <- "annotation"
    x
}

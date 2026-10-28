# mwtab2mztabm

Convert [Metabolomics Workbench](https://www.metabolomicsworkbench.org)
mwTab files to [mzTab-M](https://github.com/HUPO-PSI/mzTab) 2.1 files.

The conversion goes through two intermediate files:

1. the mwTab JSON file is downloaded from the Metabolomics Workbench,
2. its **sample metadata** is exported to a **TOML** file,
3. its **metabolite-wise quantification table** is exported to a **Parquet**
   file,
4. the mzTab-M **MTD** section is generated from the TOML file and the
   **SML** section from the Parquet file using
   [RmzTabM](https://github.com/RforMassSpectrometry/RmzTabM), and both are
   assembled into the final mzTab-M file (validated with the jmzTab-M
   validator).

## Installation

```r
install.packages("BiocManager")
BiocManager::install("integMET/mwtab2mztabm")
```

[RmzTabM](https://github.com/RforMassSpectrometry/RmzTabM) is installed from
GitHub if it is not yet available from Bioconductor.

## Quick start

```r
library(mwtab2mztabm)

## A study ID, analysis ID or the URL of the study page can be used
url <- paste0("https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?",
              "Mode=Study&StudyID=ST003934&StudyType=MS&ResultType=1")
res <- mwtabToMzTabM(url, outdir = "ST003934")
res
#>   analysis_id                          mwtab                            toml
#> 1    AN006459 ST003934/ST003934_mwtab.json ST003934/ST003934_AN006459.toml
#>                           parquet                           mztab
#> 1 ST003934/ST003934_AN006459.parquet ST003934/ST003934_AN006459.mztab
```

The study page URL returns HTML, so the study ID is extracted from it and the
mwTab JSON is downloaded from the Metabolomics Workbench REST service
(`https://www.metabolomicsworkbench.org/rest/study/study_id/ST003934/mwtab/json`).
A study can contain several analyses; one set of TOML, Parquet and mzTab-M
files is created for each.

## Step by step

```r
json <- downloadMwTab("ST003934", destdir = "ST003934")

## 1. sample metadata -> TOML
mwtabToToml(json, "ST003934/AN006459.toml", analysis_id = "AN006459")

## 2. metabolite quantification table -> Parquet
mwtabToParquet(json, "ST003934/AN006459.parquet", analysis_id = "AN006459")

## 3. TOML -> MTD section, Parquet -> SML section
tab <- readMwTabParquet("ST003934/AN006459.parquet")
mtd <- mtdFromToml("ST003934/AN006459.toml",
                   assays = attr(tab, "sample_columns"),
                   databases = smlDatabases(tab))
sml <- smlFromParquet(tab, mtd = mtd)

## 4. assemble and write the mzTab-M file (steps 3 and 4 in one call)
assembleMzTabM("ST003934/AN006459.toml", "ST003934/AN006459.parquet",
               file = "ST003934/AN006459.mztab")
```

## Intermediate files

### TOML (sample metadata)

All mwTab sections except the metabolite data are exported with their
original names and keys. The samples are an array of tables:

```toml
[mwtab2mztabm]
format_version = 1
factor_order = ["batch", "Sample source", "SRM 1950", ...]
additional_data_order = ["RAW_FILE_NAME(sample_name)", "injection_order", ...]

["METABOLOMICS WORKBENCH"]
STUDY_ID = "ST003934"
ANALYSIS_ID = "AN006459"

[STUDY]
STUDY_TITLE = "Benchmarking Data for Untargeted Mass Spectrometry-Based ..."

[MS]
INSTRUMENT_NAME = "Thermo ISQ"
ION_MODE = "POSITIVE"

[MS_METABOLITE_DATA]
Units = "not applicable"

[[SUBJECT_SAMPLE_FACTORS]]
"Subject ID" = "gcq_db5ms_pos_1_001_49.cdf"
"Sample ID" = "gcq_db5ms_pos_1_001_49.cdf"

[SUBJECT_SAMPLE_FACTORS.Factors]
batch = "1"
"Sample source" = "QC"

[SUBJECT_SAMPLE_FACTORS."Additional sample data"]
"RAW_FILE_NAME(sample_name)" = "gcq_db5ms_pos_1_001_49.cdf"
injection_order = "1"
```

TOML tables are unordered, so the original order of the factors is kept in
`[mwtab2mztabm]`.

### Parquet (quantification table)

One row per metabolite with the columns `Metabolite`, the metabolite
annotations of the mwTab `Metabolites` block (e.g. `kegg_id`, `HMDB`,
`quantified m/z`) and one `double` column per sample. The sample and
annotation column names, the units and the study/analysis IDs are stored as
key-value metadata of the file (see `?mwtabToParquet`).

## Mapping to mzTab-M

| mzTab-M | mwTab |
|---|---|
| `mzTab-ID` | `<STUDY_ID>_<ANALYSIS_ID>` |
| `title`, `description` | `STUDY_TITLE`, `STUDY_SUMMARY` |
| `contact[1]` | name, institute and email of `STUDY` |
| `uri[1]` | Metabolomics Workbench study page |
| `instrument[1]` | `MS`: `INSTRUMENT_NAME`, `MS_TYPE` (ionization), `INSTRUMENT_TYPE` (analyzer) |
| `sample[n]` | `Sample ID`, species from `SUBJECT`, `SAMPLE_TYPE` as tissue, subject ID and additional sample data as custom parameters |
| `ms_run[n]` | raw file name from the additional sample data; scan polarity from `ION_MODE` |
| `assay[n]` | one per sample column of the quantification table |
| `study_variable_group[n]`, `study_variable[n]` | sample factors and their values |
| `small_molecule-quantification_unit` | `Units` of `MS_METABOLITE_DATA` |
| SML `database_identifier`, `uri` | HMDB, KEGG or PubChem annotation column |
| SML `chemical_name` | metabolite name (only with a database identifier, see below) |
| SML `adduct_ions` | `Adduct` annotation column (e.g. `[M+H]+` becomes `[M+H]1+`) |
| SML `abundance_assay[n]` | quantification table |
| SML `abundance_study_variable[n]` | mean and coefficient of variation per study variable |
| SML `opt_global_*` | metabolite name and all annotation columns |

## Notes

- mwTab does not report the identification confidence of metabolites, but
  the mzTab-M validator requires the SML `reliability` column. It defaults
  to `"2"` (putatively annotated compound) and can be changed with
  parameter `reliability`.
- The validator requires the number of chemical names to match the number of
  database identifiers. Metabolites without an HMDB, KEGG or PubChem
  identifier therefore have `chemical_name` `null`; their mwTab name is
  available in `opt_global_metabolite_name`.
- The validator does not accept `null` as MS run location. If no raw file
  name is available, the sample ID is used.
- `RmzTabM::mzTabMValidator()` (and thus `writeMzTabM(validate = TRUE)`)
  reports files without any validation message as failed. `assembleMzTabM()`
  handles that case itself.
- Only mass spectrometry analyses with a `MS_METABOLITE_DATA` section are
  supported. Analyses whose results are in a separate file are not.

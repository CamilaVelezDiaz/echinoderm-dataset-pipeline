# Echinoderm-dataset-pipeline

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23639651)](https://doi.org/10.5281/zenodo.23639651) -->

R workflow for the integrated echinoderm occurrence dataset for northeast Australia and the adjacent Coral Sea. Produces a curated, quality-controlled, Darwin Core Archive of 43,470 echinoderm occurrence records from 13 source datasets, deposited on Zenodo and described in a data paper in Biodiversity Data Journal.

## Contents 

| File | Purpose |
|------|---------|
| `01_echinoderm_dataset_integration.R` | Loads 13 source datasets, filters to Echinodermata within the study bounding box, deduplicates records across sources, and writes long-format (Layer A) and wide-format (Layer B) tables. |
| `02_echinoderm_post_processing.R` | Resolves cross-source conflicts in coordinates, collection year, basis of record, depth, and taxonomic name (against WoRMS); assigns depth zones; generates the final `echino_wide.csv` and Darwin Core Archive. |
| `LICENSE.md` | MIT license for the code. |
| `CITATION.cff` | Machine-readable citation. |
| `sessionInfo_01.txt` | R version and package versions used to produce 01_echinoderm_dataset_integration.R v1.0.0. |
| `sessionInfo_02.txt` | R version and package versions used to produce 02_echinoderm_post_processing.R v1.0.0. |
| `.zenodo.json` | Zenodo deposit metadata. |

## Requirements

- **R version**: 4.5.1 or higher.
- **Key packages to run 01/02**: `dplyr`, `tibble`, `stringr`, `readr`, `purrr`, `tidyr`, `worrms` (WoRMS API), `sf`, `readxl`, `lubridate` and `rnaturalearth` (coordinate/land QC flagging in Section 1d of 02_echinoderm_post_processing.R - required to reproduce coord_land_qc_flag). The OBIS source is a static snapshot, generated once via `robis::occurrence()` prior to this pipeline, not regenerated at runtime, since OBIS is a continuously-updated live database and a fixed snapshot is needed for a reproducible result. robis is therefore not required to run 01_echinoderm_dataset_integration.R or 02_echinoderm_post_processing.R directly.
- Full package versions are recorded in `sessionInfo_01.txt / sessionInfo_02.txt`.

## Data sources 

The pipeline expects 13 raw source files that must be downloaded separately. Sources include:

| Number | File name |
|------|---------|
| `1` | GBIF Queensland Museum/occurrence.txt |
| `2` | GIBF - Echinodermata data_Au/occurrence.txt |
| `3` | ALA - Echinodermata dat/records-2026-04-10.csv |
| `4` | records Atlas - Other invertebrates collection/records-2026-03-13.tsv |
| `5` | Record Atlas GBRBD/records-2026-03-13.tsv |
| `6` | Csiro- GBRSBD/occurrence.txt |
| `7` | Csiro-QLDM other invertebrates - marine records/occurrence.txt |
| `8` | iDigBio data/occurrence.csv |
| `9` | OZCAM ALA data/records-2026-04-03.csv |
| `10` | OBIS/obis_echi_snapshot_2026-05-25.csv |
| `11` | Cidaris - Michela/MTQ echinoderms - registered - June 2021 final.xlsx |
| `12` | MQT - Stefano/CIDARIS_Echinoderms.xlsx |
| `13` | AM ECHINO DEPTH/26-015 Velez.csv |

The direct CMS exports supplied by QMT (11 & 12) and AM (13) were used with written permission from the depositing institutions. The original export files themselves are not redistributed through this repository; the derived occurrence records are included in the deposited dataset with permission of both institutions, acknowledged above.

## How to run

1. **Clone the repository**

    ```bash
    git clone https://github.com/CamilaVelezDiaz/echinoderm-dataset-pipeline
    cd echinoderm-dataset-pipeline
    ```

2. **Download the source datasets** into a `data-raw/` folder (see Table 3 of the *Biodiversity Data Journal* descriptor for URLs and dates).

3. **Run the pipeline**

    ```bash
    Rscript "01_echinoderm_dataset_integration.R"
    Rscript "02_echinoderm_post_processing.R"
    ```


## Main Outputs

    echino_wide.csv — the primary wide-format table (43,470 records)
    echino_wide_depth.csv — depth-bearing subset (24,684 records)
    dwca_public_release/ — Darwin Core Archive (occurrence.txt, meta.xml, eml.xml, measurementorfact.txt)

## Reproducibility notes

This pipeline was run under `LC_COLLATE=Spanish_Mexico.utf8` (see `sessionInfo_01.txt / sessionInfo_02.txt` for the full session record). This does not affect reproducibility: the only step in the pipeline where locale-dependent text sorting could influence the final result (the record_key tie-break in Section 0 of 02_echinoderm_post_processing.R`) uses locale-independent, byte-order comparison (`.locale = "C"`) rather than the session's own locale. This was verified empirically before applying the fix - zero groups in the current dataset have a genuine tie requiring this tie-break to activate - so the fix is a safeguard for future reruns with additional or updated source data, not a correction to an active problem.

Every other locale-sensitive operation in the pipeline (numeric sorting for depth/coordinate resolution, print-only diagnostic tables) either doesn't depend on text comparison at all, or is used purely for on-screen diagnostics rather than to determine which value is written to the final dataset.

Running this pipeline under any other locale (e.g. `English_US.utf8`) will produce an identical `echino_wide.csv` and Darwin Core Archive.


## Citation

If you use this code, please cite:

    Velez Diaz, M.C., Birtles, A. & Watson, S.-A. (2026). echinoderm-dataset-pipeline: R workflow for the integrated echinoderm occurrence dataset for northeast Australia and the adjacent Coral Sea (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23639651

Companion resources:

    Dataset (Darwin Core Archive): https://doi.org/10.5281/zenodo.23639651

## Licence

This code is released under the MIT License (see LICENSE).

The published integrated dataset records carry per-record licences (see license field); the authors’ compilation and derived fields are CC BY 4.0. Original CMS export files supplied by the Australian Museum and Queensland Museum Tropics remain subject to the institutions' own access conditions and are not redistributed; the derived occurrence records are included in the deposited dataset with permission of both institutions.

## Authors

    Maria Camila Velez Diaz (ORCID: 0000-0003-4180-1077) — College of Science and Engineering, James Cook University
    Alastair Birtles — College of Science and Engineering, James Cook University
    Sue-Ann Watson (ORCID: 0000-0002-9818-7429) — College of Science and Engineering, James Cook University

## Contact 
    Maria Camila Velez Diaz - mariacamila.velezdiaz@my.jcu.edu.au

## Acknowledgements
We thank the Australian Museum (Claire Rowe), Queensland Museum Tropics (Stefano Borghi) and James Cook University (Michela Mitchell) for supplying direct CMS exports of their echinoderm holdings.

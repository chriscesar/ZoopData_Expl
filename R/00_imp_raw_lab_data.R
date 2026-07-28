# 00_imp_raw_lab_data.R ----
## import lab returns for zoop data and save as an R data object

#### load packages ####
ld_pkgs <- c("tidyverse","seas","tictoc","readxl")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

tictoc::tic.clearlog()

## set universals ####
tictoc::tic("Set universals");print("set universals")
source("R/set_meta.R")
# creates object `hldfol` indicating location of raw lab returns

tictoc::toc(log=TRUE)

# find all xlsx files ----
tictoc::tic("Find all xlsx files")
xlsx_files <- list.files(
  path = hldfol,
  pattern = "\\.xlsx$",
  full.names = TRUE
  )
tictoc::toc(log=TRUE)

tictoc::tic("Extract data")
results <- purrr::map_dfr(
  xlsx_files,
  function(file) {
    
    sheets <- readxl::excel_sheets(file)
    
    purrr::map_dfr(
      sheets,
      function(sheet) {
        
        # Read metadata
        vals <- readxl::read_excel(
          path = file,
          sheet = sheet,
          range = "D2:D12",
          col_names = FALSE,
          col_types = "text"
        ) %>%
          dplyr::pull(1)
        
        spp <- readxl::read_excel(
          path = file,
          sheet = sheet,
          skip = 14,
          col_names = FALSE,
          col_types = "text"
        ) %>%
          dplyr::select(2:5)
        
        names(spp) <- c(
          "taxon",
          "aphia_id",
          "abundance_raw",
          "abundance_m3"
        )
        
        spp <- spp %>%
          dplyr::mutate(
          taxon = stringr::str_trim(taxon),
          aphia_id = stringr::str_trim(aphia_id)
          ) %>%
          dplyr::filter(!is.na(abundance_raw))
        
        # Repeat metadata for every species row
        spp %>%
          dplyr::mutate(
            filename     = basename(file),
            sheet        = sheet,
            pot_number   = vals[1],
            site_date    = vals[2],
            date = lubridate::dmy(
              stringr::str_extract(vals[2], "^\\d{2}/\\d{2}/\\d{4}")
            ),
            site = stringr::str_trim(
              stringr::str_remove(
                vals[2],
                "^\\d{2}/\\d{2}/\\d{4}"
              )),
            biosys_code  = vals[3],
            wims_code    = vals[4],
            analyst_init = vals[5],
            smp_depth_m  = vals[6],
            net_vol_m3   = vals[7],
            time         = vals[8],
            is_rep       = vals[9],
            any_other    = vals[10],
            prn          = vals[10],
            smp_comment  = vals[11],
            .before = 1
          )
      }
    )
  }
)
tictoc::toc(log=TRUE)

tictoc::tic("Write data")
write.csv(results, file = paste0(hldfol,"results.csv"),row.names = FALSE)
tictoc::toc(log=TRUE)

unlist(tictoc::tic.log())

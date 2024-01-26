library(targets)

# Source the functions and targets definitions
source("R/functions.R")

tar_option_set(packages = c("tidyverse", "magrittr", "xlsx", "readxl"))

user_excel_path <- getOption("smFishPlateDesigner_excel_path", default = NULL)


# Define the list of targets
list(
  tar_target(
    pcr_plate_layouts,
    getPCR(user_excel_path)
  ),
  tar_target(
    pcr_plate_files,
    {
      path <- "output/pcr_plate_layouts.xlsx"
      # Initialize the Excel file by writing the first plate
      write.xlsx(pcr_plate_layouts[[1]], file = path, sheetName = "Plate_1")

      # Append remaining plates to the Excel file
      if (length(pcr_plate_layouts) > 1) {
        walk2(pcr_plate_layouts[-1], seq(2, length(pcr_plate_layouts)), ~write.xlsx(.x, file = path, sheetName = paste0("Plate_", .y), append = TRUE))
      }
      path
    },
    format = "file"
  ),
  tar_target(
    pcr2_plate_layouts,
    getPCR2(pcr_plate_layouts)
  ),
  tar_target(
    export_pcr2,
    {
      path <- "output/pcr2_plate_layouts.xlsx"
      # Write each plate from pcr2_plate_layouts to a separate sheet in the Excel file
      write.xlsx(pcr2_plate_layouts[[1]], file = path, sheetName = "Plate_1")
      if (length(pcr2_plate_layouts) > 1) {
        purrr::walk2(pcr2_plate_layouts[-1], seq(2, length(pcr2_plate_layouts)), ~write.xlsx(.x, file = path, sheetName = paste0("Plate_", .y), append = TRUE))
      }
      path
    },
    format = "file"
  ),
  tar_target(
    dosage_tiv_data,
    processDosageTIV(user_excel_path)
  ),
  tar_target(
    export_dosage_tiv,
    {
      path <- "output/dosageTIV.xlsx"
      write.xlsx(dosage_tiv_data[[1]], file = path, sheetName = "Plate_1")
      if (length(dosage_tiv_data) > 1) {
        purrr::walk2(dosage_tiv_data[-1], seq(2, length(dosage_tiv_data)), ~write.xlsx(.x, file = path, sheetName = paste0("Plate_", .y), append = TRUE))
      }
      path
    },
    format = "file"
  ),
  tar_target(
    process_fish_data,
    processFishData(user_excel_path)
  ),
  tar_target(
    export_fish_data,
    {
      path <- "output/fish.xlsx"
      write.xlsx(process_fish_data[[1]], file = path, sheetName = "Plate_1")
      if (length(process_fish_data) > 1) {
        purrr::walk2(process_fish_data[-1], seq(2, length(process_fish_data)), ~write.xlsx(.x, file = path, sheetName = paste0("Plate_", .y), append = TRUE))
      }
      path
    },
    format = "file"
  ),
  tar_target(
    process_fish_data_without_primers,
    processFishDataWithoutPrimers(user_excel_path)
  ),
  tar_target(
    export_fish_data_without_primers,
    {
      path <- "output/fishWithoutPrimers.xlsx"
      write.xlsx(process_fish_data_without_primers[[1]], file = path, sheetName = "Plate_1")
      if (length(process_fish_data_without_primers) > 1) {
        purrr::walk2(process_fish_data_without_primers[-1], seq(2, length(process_fish_data_without_primers)), ~write.xlsx(.x, file = path, sheetName = paste0("Plate_", .y), append = TRUE))
      }
      path
    },
    format = "file"
  )
)

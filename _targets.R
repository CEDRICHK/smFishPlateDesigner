library(targets)

# Set Excel path before running the workflow
options(smFishPlateDesigner_excel_path = "./data/Probes_final_Oligopool11_with_sets.xlsx")

# Source the functions and targets definitions
# source("R/functions.R")

tar_option_set(packages = c("tidyverse", "magrittr", "openxlsx", "readxl", "smFishPlateDesigner"))

user_excel_path <- getOption("smFishPlateDesigner_excel_path", default = NULL)


# Define the list of targets
list(
  tar_target(
    pcr_plate_layouts,
    getPCR(user_excel_path)
  ),
  tar_target(
    pcr_plate_files,
    write_plate_workbook(pcr_plate_layouts, "output/pcr_plate_layouts.xlsx"),
    format = "file"
  ),
  tar_target(
    pcr2_plate_layouts,
    getPCR2(pcr_plate_layouts)
  ),
  tar_target(
    export_pcr2,
    write_plate_workbook(pcr2_plate_layouts, "output/pcr2_plate_layouts.xlsx"),
    format = "file"
  ),
  tar_target(
    dosage_tiv_data,
    processDosageTIV(user_excel_path)
  ),
  tar_target(
    export_dosage_tiv,
    write_plate_workbook(dosage_tiv_data, "output/dosageTIV.xlsx"),
    format = "file"
  ),
  tar_target(
    process_fish_data,
    processFishData(user_excel_path)
  ),
  tar_target(
    export_fish_data,
    write_plate_workbook(process_fish_data, "output/fish.xlsx"),
    format = "file"
  ),
  tar_target(
    process_fish_data_without_primers,
    processFishDataWithoutPrimers(user_excel_path)
  ),
  tar_target(
    export_fish_data_without_primers,
    write_plate_workbook(process_fish_data_without_primers, "output/fishWithoutPrimers.xlsx"),
    format = "file"
  )
)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----run-targets, eval=FALSE--------------------------------------------------
# # Before running the `targets` workflow, set the path to your Excel file:
# options(smFishPlateDesigner_excel_path = "path/to/your/excel_file.xlsx")
# 
# library(targets)
# tar_make()
# 
# # else No action required for the file path - uses the default embedded file
# library(targets)
# tar_make()

## ----manual-export, eval=FALSE------------------------------------------------
# plates <- getPCR("path/to/your/excel_file.xlsx")
# write_plate_workbook(plates, "output/pcr_plate_layouts.xlsx")


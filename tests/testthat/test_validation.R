test_that("validate_required_columns silently accepts complete data", {
  data <- data.frame(
    GeneName.y = "GENE1",
    BC1ID = "BC1",
    BC1PN = "PN1",
    BC1WP = "WP1",
    BC2ID = "BC2",
    BC2PN = "PN2",
    BC2WP = "WP2",
    stringsAsFactors = FALSE
  )

  expect_invisible(
    smFishPlateDesigner:::validate_required_columns(
      data,
      smFishPlateDesigner:::barcode_columns
    )
  )
})

test_that("validate_required_columns flags missing columns", {
  data <- data.frame(
    GeneName.y = "GENE1",
    BC1ID = "BC1",
    stringsAsFactors = FALSE
  )

  expect_error(
    smFishPlateDesigner:::validate_required_columns(
      data,
      smFishPlateDesigner:::barcode_columns
    ),
    "Missing required columns"
  )
})

test_that("validate_barcode_values rejects empty or NA barcodes", {
  expect_error(
    smFishPlateDesigner:::validate_barcode_values(c("A_B", NA), "unit test"),
    "Detected missing barcode values"
  )

  expect_error(
    smFishPlateDesigner:::validate_barcode_values(c("A_B", ""), "unit test"),
    "Detected missing barcode values"
  )
})

test_that("validate_barcode_values warns about duplicates", {
  out <- NULL
  expect_warning(
    out <- smFishPlateDesigner:::validate_barcode_values(
      c("A_B", "A_B", "C_D"),
      "unit test"
    ),
    "Detected duplicate barcodes"
  )
  expect_identical(out, c("A_B", "A_B", "C_D"))
})

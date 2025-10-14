make_fake_dataset <- function(n) {
  data.frame(
    GeneName.y = sprintf("GENE%03d", seq_len(n)),
    BC1ID = sprintf("BC1%03d", seq_len(n)),
    BC1PN = sprintf("PN1%03d", seq_len(n)),
    BC1WP = sprintf("WP1%03d", seq_len(n)),
    BC2ID = sprintf("BC2%03d", seq_len(n)),
    BC2PN = sprintf("PN2%03d", seq_len(n)),
    BC2WP = sprintf("WP2%03d", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

with_mocked_excel <- function(fake_data, code) {
  tmp <- tempfile(fileext = ".xlsx")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  ns <- asNamespace("smFishPlateDesigner")
  original <- get("read_excel_safe", envir = ns)
  unlockBinding("read_excel_safe", ns)
  assign("read_excel_safe", function(path, ...) fake_data, envir = ns)
  on.exit({
    assign("read_excel_safe", original, envir = ns)
    lockBinding("read_excel_safe", ns)
  }, add = TRUE)

  code(tmp)
}

test_that("getPCR generates 96-well plate layouts", {
  fake_data <- make_fake_dataset(100)

  result <- with_mocked_excel(fake_data, function(path) getPCR(path))

  expect_length(result, ceiling(100 / 96))
  expect_equal(dim(result[[1]]), c(8, 12))
  expect_identical(rownames(result[[1]]), LETTERS[1:8])
  expect_identical(colnames(result[[1]]), as.character(1:12))
})

test_that("processDosageTIV preserves ladder columns and summaries", {
  fake_data <- make_fake_dataset(120)

  result <- with_mocked_excel(fake_data, function(path) processDosageTIV(path))

  n_plate <- ceiling(120 / 96)
  expect_true("mol96" %in% names(result))
  expect_length(result$mol96, n_plate)

  plates <- result[seq_len(n_plate)]
  expect_true(all(vapply(plates, function(plate) all(plate[["12"]] == "LADDER"), logical(1))))

  bis <- result[(length(result) - n_plate + 1L):length(result)]
  expect_identical(bis[[1]][1, 8], "C-1")
  expect_true(all(vapply(bis, function(plate) all(plate[["12"]] == "LADDER"), logical(1))))
})

test_that("processFishData assembles annotated plates with controls", {
  fake_data <- make_fake_dataset(96)

  result <- with_mocked_excel(fake_data, function(path) processFishData(path))

  expect_true("mol96" %in% names(result))
  plates <- result[names(result) != "mol96"]
  expect_equal(dim(plates[[1]]), c(8, 12))
  expect_identical(plates[[1]]["F", "9"], "C-FLAP")
  expect_identical(plates[[1]]["F", "10"], "C-")
  expect_identical(plates[[1]]["G", "2"], "KIF1C")
  expect_identical(plates[[1]]["G", "3"], "DYNC1H1")
  expect_length(result$mol96, floor(length(plates) / 2))
})

test_that("processFishDataWithoutPrimers mirrors annotated layout", {
  fake_data <- make_fake_dataset(96)

  result <- with_mocked_excel(fake_data, function(path) processFishDataWithoutPrimers(path))

  expect_true("mol96" %in% names(result))
  plates <- result[names(result) != "mol96"]
  expect_equal(dim(plates[[1]]), c(8, 12))
  expect_identical(plates[[1]]["F", "9"], "C-FLAP")
  expect_identical(plates[[1]]["F", "10"], "C-")
  expect_identical(plates[[1]]["G", "2"], "KIF1C")
  expect_identical(plates[[1]]["G", "3"], "DYNC1H1")
  expect_length(result$mol96, floor(length(plates) / 2))
})

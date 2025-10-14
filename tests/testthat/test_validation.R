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

test_that("prepare_barcode_data builds unique, validated barcodes", {
  raw <- data.frame(
    GeneName.y = c("GENE1", "GENE2"),
    BC1ID = c("B1", "B2"),
    BC1PN = c("P1", "P2"),
    BC1WP = c("W1", "W2"),
    BC2ID = c("B3", "B4"),
    BC2PN = c("P3", "P4"),
    BC2WP = c("W3", "W4"),
    stringsAsFactors = FALSE
  )

  result <- smFishPlateDesigner:::prepare_barcode_data(raw)
  expect_length(result, 2)
  expect_true(all(grepl("GENE", result)))

  # Drop primer segments reduces to first element only
  result_primer <- smFishPlateDesigner:::prepare_barcode_data(
    raw,
    drop_primer_segment = TRUE
  )
  expect_identical(result_primer, c("GENE1", "GENE2"))
})

test_that("chunk_barcodes splits and pads as expected", {
  chunks <- smFishPlateDesigner:::chunk_barcodes(letters[1:5], wells_per_plate = 3, fill = "x")
  expect_length(chunks, 2)
  expect_identical(chunks[[1]], c("a", "b", "c"))
  expect_identical(chunks[[2]], c("d", "e", "x"))
})

test_that("build_plate_matrix creates plate-shaped data frame", {
  plate <- smFishPlateDesigner:::build_plate_matrix(
    chunk = LETTERS[1:6],
    nrow = 2,
    ncol = 3,
    row_labels = c("R1", "R2"),
    col_labels = c("C1", "C2", "C3")
  )

  expect_s3_class(plate, "data.frame")
  expect_identical(rownames(plate), c("R1", "R2"))
  expect_identical(colnames(plate), c("C1", "C2", "C3"))
})

test_that("build_plate_layouts wraps chunking and plate creation", {
  layout <- smFishPlateDesigner:::build_plate_layouts(
    barcodes = LETTERS[1:5],
    wells_per_plate = 4,
    nrow = 2,
    ncol = 2,
    row_labels = c("A", "B"),
    col_labels = c("1", "2"),
    fill = "-"
  )

  expect_length(layout$plates, 2)
  expect_identical(rownames(layout$plates[[1]]), c("A", "B"))
  expect_identical(colnames(layout$plates[[1]]), c("1", "2"))
  expect_identical(layout$plates[[2]][2, 2], "-")
  expect_length(layout$chunks, 2)
})

test_that("annotate_plate_set applies decorators and post-processors", {
  cfg <- list(
    wells_per_plate = 4,
    nrow = 2,
    ncol = 2,
    row_labels = c("A", "B"),
    col_labels = c("1", "2"),
    fill = "-",
    decorate = function(plate, idx) {
      plate[1, 1] <- paste0("Plate", idx)
      plate
    },
    postprocess = function(plates) {
      list(
        plates = plates,
        first_well = plates[[1]][1, 1]
      )
    }
  )

  annotated <- smFishPlateDesigner:::annotate_plate_set(
    barcodes = LETTERS[1:5],
    layout_config = cfg
  )

  expect_named(annotated, c("plates", "first_well"))
  expect_identical(annotated$first_well, "Plate1")
  expect_identical(rownames(annotated$plates[[1]]), c("A", "B"))
})

test_that("compose_plate_decorator injects fixed controls", {
  decorator <- smFishPlateDesigner:::compose_plate_decorator(
    controls = list(
      fixed = list("A,1" = "CTRL"),
      callbacks = list(function(plate, idx) {
        plate[1, 2] <- paste0("IDX", idx)
        plate
      })
    )
  )

  plate <- data.frame(
    `1` = c("a", "c"),
    `2` = c("b", "d"),
    row.names = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  decorated <- decorator(plate, 5)
  expect_identical(decorated["A", "1"], "CTRL")
  expect_identical(decorated["A", "2"], "IDX5")
})

test_that("dosage_tiv_layout_config postprocess mirrors legacy structure", {
  cfg <- smFishPlateDesigner:::dosage_tiv_layout_config()
  barcodes <- sprintf("BC%03d", seq_len(cfg$wells_per_plate))

  layout <- smFishPlateDesigner:::build_plate_layouts(
    barcodes = barcodes,
    wells_per_plate = cfg$wells_per_plate,
    nrow = cfg$nrow,
    ncol = cfg$ncol,
    row_labels = cfg$row_labels,
    col_labels = cfg$col_labels
  )

  result <- cfg$postprocess(
    plates = layout$plates,
    chunks = layout$chunks,
    config = cfg
  )

  expect_length(result, 3)
  expect_true(all(result[[1]][[cfg$col_labels[cfg$ncol]]] == "LADDER"))
  expect_identical(result[[2]], sprintf("BC%03d", cfg$wells_per_plate))
  expect_true(all(result[[3]][[cfg$col_labels[cfg$ncol]]] == "LADDER"))
})

test_that("fish_layout_config regenerates annotated plates", {
  cfg <- smFishPlateDesigner:::fish_layout_config()
  barcodes <- sprintf("FB%03d", seq_len(96))

  layout <- smFishPlateDesigner:::build_plate_layouts(
    barcodes = barcodes,
    wells_per_plate = cfg$wells_per_plate,
    nrow = cfg$nrow,
    ncol = cfg$ncol,
    row_labels = cfg$row_labels,
    col_labels = cfg$col_labels,
    decorate = cfg$decorate
  )

  result <- cfg$postprocess(
    plates = layout$plates,
    chunks = layout$chunks,
    config = cfg
  )

  expect_length(result, length(layout$plates) + 1)
  expect_identical(rownames(result[[1]]), LETTERS[1:8])
  expect_identical(colnames(result[[1]]), as.character(seq_len(12)))
  expect_identical(result[[1]][7, 2], "KIF1C")
  expect_identical(result[[1]][7, 3], "DYNC1H1")
  expect_identical(result[[1]][6, 9], "C-FLAP")
  expect_identical(result[[1]][6, 10], "C-")
  expect_length(result[[length(result)]], floor(length(layout$plates) / 2))
})

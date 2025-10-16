#' Validate or resolve the path to an Excel workbook
#'
#' Normalises the optional user-supplied `file_path` by falling back to the
#' package demo workbook and checks the file exists.
#'
#' @param file_path Path supplied by the caller. May be `NULL`.
#' @param default Path to the bundled demo workbook.
#' @return Character scalar with the resolved, readable path.
#' @keywords internal
validate_excel_path <- function(file_path,
                                default = system.file("extdata", "ResultWithComplementSeq_final_file.xlsx",
                                                      package = "smFishPlateDesigner")) {
  if (is.null(file_path)) {
    file_path <- default
  }
  if (is.null(file_path) || file_path == "") {
    stop("No Excel file was provided and no package default could be resolved.")
  }
  if (!fs::file_exists(file_path)) {
    stop("Unable to find the Excel file: ", file_path)
  }
  file_path
}

#' Safely read an Excel workbook with an informative error
#'
#' Executes `readxl::read_excel()` inside a `tryCatch()` to surface precise
#' diagnostics when the import fails (missing sheet, corrupt file, etc.).
#'
#' @param path Path to the workbook to import.
#' @param ... Additional arguments passed to `readxl::read_excel()`.
#' @return A tibble containing the workbook contents.
#' @keywords internal
read_excel_safe <- function(path, ...) {
  tryCatch(
    readxl::read_excel(path = path, ...),
    error = function(err) {
      stop(
        "Failed to read Excel file '", path, "': ",
        conditionMessage(err),
        call. = FALSE
      )
    }
  )
}

#' Ensure an Excel worksheet exposes the expected columns
#'
#' Validates that `data` contains the complete set of `required` column names
#' before downstream processing continues.
#'
#' @param data Data frame or tibble read from the workbook.
#' @param required Character vector of columns that must be present.
#' @return Invisibly returns `data` when validation succeeds.
#' @keywords internal
validate_required_columns <- function(data, required) {
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(data)
}

#' Default set of columns used to construct barcodes
#'
#' Shared helper constant listing the Excel columns that are concatenated to
#' form the unique barcode identifiers across the package.
#' @keywords internal
barcode_columns <- c(
  "GeneName.y", "BC1ID", "BC1PN", "BC1WP",
  "BC2ID", "BC2PN", "BC2WP"
)

#' Flag missing or duplicated barcode identifiers
#'
#' Ensures barcode vectors do not contain missing values and surfaces
#' duplicates so the caller can inspect the input data set.
#'
#' @param barcodes Character vector of barcode identifiers.
#' @param context Short label describing the calling function for diagnostics.
#' @return Invisibly returns `barcodes` when validation succeeds.
#' @keywords internal
validate_barcode_values <- function(barcodes, context) {
  if (anyNA(barcodes) || any(barcodes == "")) {
    stop("Detected missing barcode values in ", context, ".", call. = FALSE)
  }

  duplicates <- unique(barcodes[duplicated(barcodes)])
  if (length(duplicates) > 0L) {
    warning(
      "Detected duplicate barcodes in ", context, ": ",
      paste(head(duplicates, 10), collapse = ", "),
      if (length(duplicates) > 10) " ...",
      call. = FALSE
    )
  }

  invisible(barcodes)
}

#' Prepare a vector of barcodes from a raw data frame
#'
#' Selects the expected columns, concatenates them with an underscore, and
#' optionally strips primer information before returning a (optionally unique)
#' character vector of barcodes.
#'
#' @param data Data frame read from the Excel workbook.
#' @param columns Character vector of columns to concatenate.
#' @param drop_primer_segment Logical; if `TRUE`, keep only the first segment of
#'   each barcode (used when primers must be removed).
#' @param unique_only Logical; deduplicate barcodes when `TRUE`.
#' @param sep Separator used when concatenating the columns.
#' @param context Label describing the caller for diagnostic messages.
#' @return Character vector of barcodes.
#' @keywords internal
prepare_barcode_data <- function(data,
                                 columns = barcode_columns,
                                 drop_primer_segment = FALSE,
                                 unique_only = TRUE,
                                 sep = "_",
                                 context = "prepare_barcode_data()") {
  selected <- dplyr::select(data, dplyr::all_of(columns))

  if (anyNA(selected)) {
    stop("Detected missing values in barcode source columns.", call. = FALSE)
  }

  barcode_tbl <- tidyr::unite(
    selected,
    col = "BC",
    dplyr::everything(),
    sep = sep,
    remove = FALSE
  )

  barcodes <- barcode_tbl$BC

  if (drop_primer_segment) {
    barcodes <- vapply(
      strsplit(barcodes, sep, fixed = TRUE),
      `[`,
      FUN.VALUE = character(1),
      1
    )
  }

  validate_barcode_values(barcodes, context)

  if (unique_only) {
    barcodes <- unique(barcodes)
  }

  barcodes
}

#' Split a barcode vector into plate-sized chunks
#'
#' Pads the barcode vector to a multiple of `wells_per_plate` and splits it into
#' a list of equal-length vectors, each representing a plate.
#'
#' @param barcodes Character vector of barcodes.
#' @param wells_per_plate Number of wells on a plate.
#' @param fill Value used to pad incomplete plates.
#' @return List of barcode vectors, one per plate.
#' @keywords internal
chunk_barcodes <- function(barcodes,
                           wells_per_plate = 96,
                           fill = NA_character_) {
  if (!length(barcodes)) {
    return(list())
  }

  if (!is.numeric(wells_per_plate) || wells_per_plate <= 0) {
    stop("`wells_per_plate` must be a positive integer.", call. = FALSE)
  }

  barcodes <- as.character(barcodes)
  n_plate <- ceiling(length(barcodes) / wells_per_plate)
  total <- n_plate * wells_per_plate
  if (total > length(barcodes)) {
    barcodes <- c(barcodes, rep(fill, total - length(barcodes)))
  }

  split(barcodes, rep(seq_len(n_plate), each = wells_per_plate))
}

#' Convert a barcode chunk into a plate data frame
#'
#' Turns a character vector into a plate-shaped data frame with optional row and
#' column labels and an optional decoration function that can alter the plate
#' content (for control wells, etc.).
#'
#' @param chunk Character vector representing the barcodes for a single plate.
#' @param nrow Number of rows of the plate.
#' @param ncol Number of columns of the plate.
#' @param row_labels Optional row names.
#' @param col_labels Optional column names.
#' @param byrow Logical; fill the matrix by row when `TRUE`.
#' @param decorate Optional function called with `(plate_df, plate_index)` to
#'   modify the plate after creation.
#' @param plate_index Index of the plate (passed to `decorate`).
#' @return Data frame shaped like the plate.
#' @keywords internal
build_plate_matrix <- function(chunk,
                               nrow,
                               ncol,
                               row_labels = LETTERS[seq_len(nrow)],
                               col_labels = seq_len(ncol),
                               byrow = TRUE,
                               decorate = NULL,
                               plate_index = 1) {
  expected <- nrow * ncol
  if (length(chunk) < expected) {
    chunk <- c(chunk, rep(NA_character_, expected - length(chunk)))
  } else if (length(chunk) > expected) {
    stop(
      "Chunk length does not match expected plate size (",
      expected,
      ").",
      call. = FALSE
    )
  }

  plate_matrix <- matrix(chunk, nrow = nrow, ncol = ncol, byrow = byrow)
  colnames(plate_matrix) <- col_labels
  rownames(plate_matrix) <- row_labels

  plate_df <- as.data.frame(
    plate_matrix,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!is.null(decorate)) {
    plate_df <- decorate(plate_df, plate_index)
  }

  plate_df
}

#' Build plate layouts from a barcode vector
#'
#' Convenience wrapper that chunks the barcode vector and converts each chunk
#' into a plate-shaped data frame, optionally decorating each plate.
#'
#' @param barcodes Character vector of barcodes.
#' @param wells_per_plate Number of wells per plate.
#' @param nrow Number of rows in the plate.
#' @param ncol Number of columns in the plate.
#' @param row_labels Optional row names.
#' @param col_labels Optional column names.
#' @param byrow Logical; fill the matrix by row when `TRUE`.
#' @param decorate Optional function applied to each plate.
#' @param fill Padding value for incomplete plates.
#' @return List of plate data frames.
#' @keywords internal
build_plate_layouts <- function(barcodes,
                                wells_per_plate = 96,
                                nrow = 8,
                                ncol = 12,
                                row_labels = LETTERS[seq_len(nrow)],
                                col_labels = as.character(seq_len(ncol)),
                                byrow = TRUE,
                                decorate = NULL,
                                fill = NA_character_) {
  chunks <- chunk_barcodes(barcodes, wells_per_plate, fill = fill)

  plates <- lapply(seq_along(chunks), function(idx) {
    build_plate_matrix(
      chunk = chunks[[idx]],
      nrow = nrow,
      ncol = ncol,
      row_labels = row_labels,
      col_labels = col_labels,
      byrow = byrow,
      decorate = decorate,
      plate_index = idx
    )
  })

  list(
    plates = plates,
    chunks = chunks
  )
}

#' Assemble and post-process plate layouts according to a configuration
#'
#' Serves as a high-level orchestrator that builds plate layouts from a barcode
#' vector and applies optional decoration/post-processing steps defined in a
#' configuration list.
#'
#' Expected elements in `layout_config`:
#' - `wells_per_plate` (numeric, default 96)
#' - `nrow`, `ncol` (numeric, default 8x12)
#' - `row_labels`, `col_labels` (character vectors)
#' - `byrow` (logical)
#' - `decorate` (function `function(plate_df, plate_index)` applied to each plate)
#' - `fill` (padding value, default `NA_character_`)
#' - `postprocess` (function `function(plates, chunks, config)` returning the final result)
#' - `controls` (list describing control placement rules)
#'
#' The optional `controls` list accepts:
#' - `fixed`: named list whose names are `"ROW,COL"` coordinates (e.g. `"H,12"`)
#'   mapping to values to set in the corresponding wells.
#' - `callbacks`: list of functions `function(plate_df, plate_index)` executed
#'   in order to apply more advanced control logic.
#'
#' @param barcodes Character vector of barcodes.
#' @param layout_config List describing how plates should be constructed.
#' @return Result of `postprocess(plates, chunks, cfg)` when supplied, otherwise the list of plates.
#' @keywords internal
annotate_plate_set <- function(barcodes, layout_config = list()) {
  defaults <- list(
    wells_per_plate = 96,
    nrow = 8,
    ncol = 12,
    row_labels = LETTERS[1:8],
    col_labels = as.character(seq_len(12)),
    byrow = TRUE,
    decorate = NULL,
    fill = NA_character_,
    postprocess = NULL,
    controls = NULL
  )

  cfg <- utils::modifyList(defaults, layout_config, keep.null = TRUE)

  decorate <- cfg$decorate
  if (!is.null(cfg$controls)) {
    decorate <- compose_plate_decorator(decorate, cfg$controls)
  }

  layout <- build_plate_layouts(
    barcodes = barcodes,
    wells_per_plate = cfg$wells_per_plate,
    nrow = cfg$nrow,
    ncol = cfg$ncol,
    row_labels = cfg$row_labels,
    col_labels = cfg$col_labels,
    byrow = cfg$byrow,
    decorate = decorate,
    fill = cfg$fill
  )

  if (is.function(cfg$postprocess)) {
    args <- list(
      plates = layout$plates,
      chunks = layout$chunks,
      config = cfg
    )
    formal_names <- names(formals(cfg$postprocess))
    if (is.null(formal_names)) {
      formal_names <- character()
    }
    if (!"..." %in% formal_names) {
      args <- args[intersect(names(args), formal_names)]
    }
    return(do.call(cfg$postprocess, args))
  }

  layout$plates
}
#' Combine control placement and decoration for a plate
#'
#' Creates a decorator that first applies fixed control values and callback
#' functions before delegating to an optional downstream decorator.
#'
#' @param decorate Optional base decorator function.
#' @param controls Optional list describing controls (see `annotate_plate_set()`).
#' @return A decorator function accepting `(plate_df, plate_index)`.
#' @keywords internal
compose_plate_decorator <- function(decorate = NULL, controls = NULL) {
  fixed_controls <- NULL
  callbacks <- NULL

  if (!is.null(controls)) {
    if (!is.null(controls$fixed)) {
      fixed_controls <- controls$fixed
    }
    if (!is.null(controls$callbacks)) {
      callbacks <- controls$callbacks
    }
  }

  function(plate, idx) {
    if (!is.null(fixed_controls)) {
      for (coord in names(fixed_controls)) {
        pieces <- strsplit(coord, ",", fixed = TRUE)[[1]]
        if (length(pieces) != 2) {
          stop("Invalid control coordinate: ", coord, call. = FALSE)
        }
        plate[pieces[1], pieces[2]] <- fixed_controls[[coord]]
      }
    }

    if (is.list(callbacks) && length(callbacks)) {
      for (fn in callbacks) {
        plate <- fn(plate, idx)
      }
    }

    if (is.function(decorate)) {
      plate <- decorate(plate, idx)
    }

    plate
  }
}

#' Layout configuration for dosage TIV plate exports
#'
#' Provides the parameters and post-processing logic required to build dosage
#' TIV plate layouts using `annotate_plate_set()`.
#' @return A list describing the layout configuration.
#' @keywords internal
dosage_tiv_layout_config <- function() {
  list(
    wells_per_plate = 96,
    nrow = 8,
    ncol = 12,
    row_labels = LETTERS[1:8],
    col_labels = as.character(seq_len(12)),
    postprocess = function(plates, chunks, config) {
      if (!length(plates)) {
        return(list())
      }

      last_col <- config$col_labels[config$ncol]

      col12_values <- lapply(
        plates,
        function(plate) as.character(plate[[last_col]])
      )

      plates_with_ladder <- lapply(plates, function(plate) {
        plate[[last_col]] <- "LADDER"
        plate
      })

      padded_col12 <- lapply(col12_values, function(vals) {
        c(vals, rep(NA_character_, config$wells_per_plate - length(vals)))
      })

      bis <- lapply(seq_along(padded_col12), function(i) {
        mat <- matrix(
          padded_col12[[i]],
          nrow = config$nrow,
          ncol = config$ncol,
          byrow = TRUE
        )
        mat[1, 8] <- paste0("C-", i)
        df <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
        colnames(df) <- config$col_labels
        rownames(df) <- config$row_labels
        df[[last_col]] <- "LADDER"
        df
      })

      mol96 <- unname(vapply(
        chunks,
        function(chunk) as.character(chunk[config$wells_per_plate]),
        character(1)
      ))

      c(plates_with_ladder, list(mol96 = mol96), bis)
    }
  )
}

#' Layout configuration for fish plate exports
#'
#' Supplies the parameters and post-processing rules used by the fish workflows
#' (`processFishData()` and `processFishDataWithoutPrimers()`).
#' @return A list describing the layout configuration.
#' @keywords internal
fish_layout_config <- function() {
  final_cols <- as.character(seq_len(12))
  final_rows <- LETTERS[1:8]

  apply_base_controls <- function(plate, plate_index) {
    plate[5, 9] <- "C-FLAP"
    plate[5, 10] <- "C-"
    plate
  }

  apply_extended_controls <- function(plate, idx, kif_col, dync_col) {
    original_f9 <- plate["F", "9"]

    plate["F", "10"] <- "C-FLAP"
    plate["F", "11"] <- "C-"
    plate["F", "12"] <- NA_character_

    plate[7, kif_col] <- "KIF1C"
    plate[7, dync_col] <- "DYNC1H1"

    if (idx %% 2 == 0) {
      plate["F", "9"] <- paste0("C-", idx / 2)
    }

    list(
      plate = plate,
      mol96_entry = original_f9
    )
  }

  list(
    wells_per_plate = 48,
    nrow = 5,
    ncol = 10,
    row_labels = LETTERS[2:6],
    col_labels = as.character(2:11),
    decorate = function(plate, idx) apply_base_controls(plate, idx),
    postprocess = function(plates, chunks, config) {
      if (!length(plates)) {
        return(list())
      }

      kif_col <- 2L
      dync_col <- 3L
      mol96 <- character()
      plates_final <- vector("list", length(plates))

      for (idx in seq_along(plates)) {
        plate <- plates[[idx]]
        plate[] <- lapply(plate, as.character)

        plate <- cbind("1" = NA_character_, plate)
        plate[["12"]] <- NA_character_
        plate <- plate[, final_cols, drop = FALSE]

        extra <- matrix(
          NA_character_,
          nrow = 3,
          ncol = length(final_cols),
          dimnames = list(c("A", "G", "H"), final_cols)
        )
        plate <- rbind(plate, extra)
        rownames(plate) <- c(LETTERS[2:6], "A", "G", "H")
        plate <- plate[order(rownames(plate)), , drop = FALSE]
        rownames(plate) <- final_rows
        plate[] <- lapply(plate, as.character)

        decorated <- apply_extended_controls(plate, idx, kif_col, dync_col)
        plate <- decorated$plate

        kif_col <- kif_col + 1L
        dync_col <- dync_col + 1L
        if (kif_col == 11L) {
          kif_col <- 2L
          dync_col <- 3L
        }

        if (idx %% 2 == 0) {
          mol96 <- c(mol96, decorated$mol96_entry)
        }

        plates_final[[idx]] <- plate
      }

      c(plates_final, list(mol96 = mol96))
    }
  )
}

#' Write a list of plates to an Excel workbook
#'
#' Serializes a list of tabular objects (plate layouts) into an Excel workbook,
#' one sheet per plate. Pure vectors are coerced into one-column data frames
#' using the sheet name (when available) as the column name.
#'
#' @param plates List of plate layouts (data frames, matrices, or vectors).
#' @param path Destination file path for the workbook.
#' @param prefix Prefix used to generate sheet names when none are provided.
#' @param sheet_namer Optional callback `function(plate, index)` returning the
#'   sheet name for each plate.
#' @return The `path` of the created workbook.
#' @export
write_plate_workbook <- function(plates,
                                 path,
                                 prefix = "Plate_",
                                 sheet_namer = NULL) {
  if (length(plates) == 0L) {
    stop("No plates supplied for export.", call. = FALSE)
  }

  dir_path <- dirname(path)
  if (!identical(dir_path, ".") && !fs::dir_exists(dir_path)) {
    fs::dir_create(dir_path, recurse = TRUE)
  }

  if (fs::file_exists(path)) {
    fs::file_delete(path)
  }

  sheet_counter <- 0L
  plate_names <- names(plates)
  wb <- openxlsx::createWorkbook()

  for (idx in seq_along(plates)) {
    plate_object <- plates[[idx]]
    if (is.null(plate_object)) {
      next
    }

    original_is_vector <- is.null(dim(plate_object)) && !is.list(plate_object)
    original_has_rownames <- !is.null(rownames(plate_object))

    if (original_is_vector) {
      plate_df <- data.frame(
        value = unname(plate_object),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      original_has_rownames <- FALSE
    } else if (is.matrix(plate_object)) {
      plate_df <- as.data.frame(plate_object, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      plate_df <- as.data.frame(plate_object, stringsAsFactors = FALSE, check.names = FALSE)
    }

    if (!is.data.frame(plate_df)) {
      next
    }

    if (original_has_rownames) {
      row_col <- rownames(plate_df)
      plate_df <- cbind(row_col, plate_df)
      colnames(plate_df)[1] <- ""
    }

    sheet_counter <- sheet_counter + 1L
    default_name <- paste0(prefix, sheet_counter)

    sheet_name <- if (!is.null(sheet_namer)) {
      sheet_namer(plate_df, idx)
    } else {
      candidate <- plate_names[idx]
      if (original_is_vector) {
        if (!is.null(candidate) && nzchar(candidate)) {
          candidate <- candidate
        } else {
          candidate <- NULL
        }
      }
      if (!is.null(candidate) && nzchar(candidate)){
        candidate
      } else {
        default_name
      }
    }

    if (is.null(sheet_name) || !nzchar(sheet_name)) {
      sheet_name <- default_name
    }

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = plate_df,
      colNames = TRUE,
      rowNames = FALSE
    )
  }

  if (sheet_counter == 0L) {
    stop("No tabular data found to export.", call. = FALSE)
  }

  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)

  path
}


#' Create Plate Layout for PCR from User-Provided Excel File
#'
#' This function processes PCR data from a user-provided Excel file and organizes the data into
#' plate layouts. It assumes each plate has 96 wells and arranges unique barcode sequences
#' into an 8x12 matrix corresponding to the layout of a standard PCR plate.
#'
#' @param file_path A string specifying the path to the Excel file containing PCR data.
#'                  The file should be in a format readable by `readxl::read_excel()`.
#' @return A list of data frames, each representing the layout of a PCR plate.
#'         Rows are labeled A-H and columns 1-12, with each cell containing a unique barcode.
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "ResultWithComplementSeq_final_file.xlsx", package = "smFishPlateDesigner")
#' plates <- getPCR(file_path)
#' }
#' @export
#'
#' @details
#' Internally the function delegates to shared helpers that validate the Excel
#' schema, create a unique set of barcodes via `prepare_barcode_data()`, then
#' assemble a list of 96-well plate layouts with `annotate_plate_set()`. Any
#' missing barcode values raise an error while duplicates generate a warning.
#' @importFrom readxl read_excel
getPCR <- function(file_path=NULL) {

  file_path <- validate_excel_path(file_path)

  # Read the Excel file at the given file path
  data <- read_excel_safe(file_path)
  validate_required_columns(data, barcode_columns)

  barcodes <- prepare_barcode_data(
    data,
    columns = barcode_columns,
    drop_primer_segment = FALSE,
    unique_only = TRUE,
    context = "getPCR() input"
  )

  annotate_plate_set(
    barcodes,
    layout_config = list(
      wells_per_plate = 96,
      nrow = 8,
      ncol = 12,
      row_labels = LETTERS[1:8],
      col_labels = as.character(seq_len(12))
    )
  )
}

#' Further Process PCR Data and Annotate Plate Matrices
#'
#' Takes the output from `getPCR` function, performs additional processing,
#' extracts gene names from the 96th position, and annotates the plate matrices with "C-".
#'
#' @param pcr_data A list of data frames representing plate layouts from `getPCR`.
#' @return A list containing the annotated plate data frames and a vector of gene names from the 96th position.
#' @examples
#' \dontrun{
#' # Example assumes `pcr_data` is available from `getPCR`.
#' getPCR2(pcr_data)
#' }
#' @export
getPCR2 <- function(pcr_data) {

  # Extract the gene name from position H12 (96th position) of each plate
  mol96 <- sapply(pcr_data, function(x) {
    as.character(x[8, 12])
  })

  # Annotate each plate with "C-" and the plate number at position H12
  l_mat_annotated <- lapply(seq_along(pcr_data), function(x) {
    plate <- pcr_data[[x]]
    plate <- plate %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    plate[8, 12] <- paste("C-", x, sep = "")
    plate
  })

  # Combine the list of annotated plates with the control list of gene names
  lctrls <- list(mol96 = mol96)
  l_annotated <- c(l_mat_annotated, lctrls)

  return(l_annotated)
}


#' Process Dosage TIV Data from User-Provided Excel File
#'
#' This function processes TIV dosage data from an Excel file specified by the user.
#' The user must provide the path to their Excel file containing the dosage data.
#'
#' @param file_path A string specifying the path to the Excel file to be processed.
#'                  The file should be in a format readable by `readxl::read_excel()`.
#' @return A data frame containing the processed dosage data.
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "ResultWithComplementSeq_final_file.xlsx", package = "smFishPlateDesigner")
#' processed_data <- processDosageTIV(file_path)
#' }
#' @export
#'
#' @importFrom readxl read_excel
#' @details
#' Uses the shared helper pipeline (`prepare_barcode_data()` +
#' `annotate_plate_set()`) to clean the Excel input and construct plate layouts.
#' Column 12 is overwritten with the `LADDER` marker, a secondary "bis" workbook
#' mirroring the legacy output is appended, and the identifiers from the H12
#' well of every plate are returned as the `mol96` element.
processDosageTIV <- function(file_path = NULL) {
  file_path <- validate_excel_path(file_path)
  data <- read_excel_safe(file_path)
  validate_required_columns(data, barcode_columns)

  barcodes <- prepare_barcode_data(
    data,
    columns = barcode_columns,
    drop_primer_segment = FALSE,
    unique_only = TRUE,
    context = "processDosageTIV() input"
  )

  annotate_plate_set(
    barcodes,
    layout_config = dosage_tiv_layout_config()
  )
}



#' Process Fish Data into Plate Layouts from User-Provided Excel File
#'
#' This function processes fish data from a user-provided Excel file, creating unique barcodes for each sample
#' and arranging them into plate layouts with specific annotations. It is designed to read the Excel file and
#' perform necessary processing to generate a structured layout for laboratory use.
#'
#' @param file_path A string specifying the path to the Excel file containing fish data.
#'                  The file should be in a format readable by `readxl::read_excel()`.
#' @return A list containing the plate layouts with annotations. Each layout will correspond
#'         to the structure of a standard PCR plate, with specific barcode allocations.
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "ResultWithComplementSeq_final_file.xlsx", package = "smFishPlateDesigner")
#' processFishData(file_path)
#' }
#' @export
#'
#' @details
#' The function reads the specified Excel file and funnels it through the shared
#' validation/build pipeline, ultimately delegating to `annotate_plate_set()`
#' with the fish-specific layout configuration. The configuration pads the base
#' 5x10 layout to a full 96-well plate, inserts control markers (`C-FLAP`,
#' `C-`, `KIF1C`, `DYNC1H1`), and labels every second plate with its `C-`
#' sample identifier in well H9. Duplicate barcodes still trigger warnings while
#' missing barcodes abort the processing.
#' @importFrom readxl read_excel
processFishData <- function(file_path=NULL) {
  file_path <- validate_excel_path(file_path)
  # Import data file and process
  data <- read_excel_safe(file_path)
  validate_required_columns(data, barcode_columns)

  barcodes <- prepare_barcode_data(
    data,
    columns = barcode_columns,
    drop_primer_segment = FALSE,
    unique_only = TRUE,
    context = "processFishData() input"
  )

  if (length(barcodes)) {
    nblock96 <- ceiling(length(barcodes) / 96)
    target_length <- nblock96 * 96
    if (target_length > length(barcodes)) {
      barcodes <- c(
        barcodes,
        rep(NA_character_, target_length - length(barcodes))
      )
    }
  }

  annotate_plate_set(
    barcodes,
    layout_config = fish_layout_config()
  )
}


#' Process Fish Data into Plate Layouts from User-Provided Excel File
#'
#' This function processes fish data from a user-provided Excel file, creating unique barcodes for each sample
#' and arranging them into plate layouts with specific annotations. The function is tailored to handle data
#' where primer information is not required or is pre-processed.
#'
#' @param file_path A string specifying the path to the Excel file containing fish data.
#'                  The file should be in a format readable by `readxl::read_excel()`.
#' @return A list containing the plate layouts with annotations. Each layout will correspond
#'         to the structure of a standard PCR plate, with specific barcode allocations.
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "ResultWithComplementSeq_final_file.xlsx", package = "smFishPlateDesigner")
#' processFishDataWithoutPrimers(file_path)
#' }
#' @export
#'
#' @details
#' The function reads the specified Excel file, validates the expected columns,
#' strips primer segments through `prepare_barcode_data(drop_primer_segment = TRUE)`,
#' and delegates plate construction to the shared fish configuration used by
#' `processFishData()`. The same control markers and `mol96` summary element are
#' produced, ensuring behaviour identical to the legacy implementation.
#' @importFrom readxl read_excel
processFishDataWithoutPrimers <- function(file_path=NULL) {
  file_path <- validate_excel_path(file_path)
  #data <- read_excel(path = args[5])
  data <- read_excel_safe(file_path)
  validate_required_columns(data, barcode_columns)

  barcodes <- prepare_barcode_data(
    data,
    columns = barcode_columns,
    drop_primer_segment = TRUE,
    unique_only = TRUE,
    context = "processFishDataWithoutPrimers() input"
  )

  if (length(barcodes)) {
    nblock96 <- ceiling(length(barcodes) / 96)
    target_length <- nblock96 * 96
    if (target_length > length(barcodes)) {
      barcodes <- c(
        barcodes,
        rep(NA_character_, target_length - length(barcodes))
      )
    }
  }

  annotate_plate_set(
    barcodes,
    layout_config = fish_layout_config()
  )
}

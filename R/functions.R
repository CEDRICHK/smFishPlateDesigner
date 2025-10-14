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
      if (length(duplicates) > 10) " …",
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
                                col_labels = seq_len(ncol),
                                byrow = TRUE,
                                decorate = NULL,
                                fill = NA_character_) {
  chunks <- chunk_barcodes(barcodes, wells_per_plate, fill = fill)

  lapply(seq_along(chunks), function(idx) {
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
#' - `postprocess` (function `function(plates)` returning the final result)
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
#' @return Result of `postprocess(plates)` when supplied, otherwise the list of plates.
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

  plates <- build_plate_layouts(
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
    return(cfg$postprocess(plates))
  }

  plates
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
#' The function reads the user-provided Excel file, verifies that all required
#' barcode columns are present, and fails fast if any are missing. It then
#' creates a unified barcode by concatenating these columns with an underscore
#' separator. If any barcode values are empty or `NA`, the function aborts with
#' a clear diagnostic, and duplicated barcodes trigger a warning. After
#' deduplicating barcodes, it calculates the number of plates needed and
#' arranges barcodes into an 8x12 matrix for each plate. The final list of
#' matrices can be used to guide the setup of PCR plates in a laboratory
#' setting.
#' @importFrom readxl read_excel
getPCR <- function(file_path=NULL) {

  file_path <- validate_excel_path(file_path)

  # Read the Excel file at the given file path
  data <- read_excel_safe(file_path)
  validate_required_columns(data, barcode_columns)

  # Selecting relevant columns for barcode creation
  selected_data <- data %>%
    dplyr::select(GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP)

  # Creating a unified barcode by concatenating selected columns with an underscore
  barcode <- selected_data %>%
    tidyr::unite(col = BC, c(GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP), sep = "_")



  validate_barcode_values(barcode$BC, "getPCR() input")

  # Removing duplicate barcodes
  unique_barcode <- dplyr::distinct(barcode, BC)

  # Calculating the number of 96-well plates needed based on the number of unique barcodes
  n_plate <- ceiling(nrow(unique_barcode) / 96)

  # Ensuring the number of barcodes fits into the 96-well plate format by adding empty wells if necessary
  barcode_vector <- c(unique_barcode$BC, rep(NA, 96 * n_plate - nrow(unique_barcode)))

  # Creating a list of matrices, each representing a PCR plate
  plates <- lapply(seq_len(n_plate), function(i) {
    plate_matrix <- matrix(
      barcode_vector[(((i - 1) * 96) + 1):(i * 96)],
      nrow = 8,
      ncol = 12,
      byrow = TRUE
    )
    # Annotating the plate matrix with row (A-H) and column (1-12) labels
    colnames(plate_matrix) <- 1:12
    rownames(plate_matrix) <- LETTERS[1:8]
    as.data.frame(plate_matrix) # Converting matrix to a data frame for better readability
  })

  return(plates)
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
    plate %<>% mutate_all(as.character)
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
#' Validates the presence of all required barcode columns before processing.
#' The function aborts if any barcode value is empty or `NA`, and it emits a
#' warning when duplicate barcodes are encountered. On success, it returns the
#' annotated plate layouts derived from the cleaned barcode set.
processDosageTIV <- function(file_path = NULL) {
  file_path <- validate_excel_path(file_path)
  data <- read_excel_safe(file_path)
  validate_required_columns(data, barcode_columns)

  # Select features and create a barcode
  data <- dplyr::select(data, GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP)
  barcode <- tidyr::unite(data, col = BC, GeneName.y, BC1ID, BC1PN,
                          BC1WP, BC2ID, BC2PN, BC2WP, sep = "_")

  validate_barcode_values(barcode$BC, "processDosageTIV() input")

  # Process barcode and calculate number of plates
  u_barcode <- dplyr::distinct(barcode, BC, .keep_all = TRUE)$BC
  nplate <- ceiling(length(u_barcode) / 96)
  u_barcode <- append(u_barcode, rep(NA, 96 * nplate - length(u_barcode)))

  # Create and annotate the matrices
  l_mat <- lapply(seq_len(nplate), function(x) {
    matrix(u_barcode[((x - 1) * 96 + 1):(x * 96)], nrow = 8, byrow = TRUE)
  })

  l_matannotated <- lapply(l_mat, function(mat) {
    df <- as.data.frame(mat)
    colnames(df) <- 1:12
    rownames(df) <- LETTERS[1:8]
    df
  })

  #get name on position 96 and replace value by "C-1"
  v <- vector()
  mol96 <- sapply(l_matannotated, function(x) {
    v <- append(v, as.character(x[8,12]))
    return(v)
  })

  #get names of the column 12 and replace value by "LADDER"
  v <- vector()
  col12 <- lapply(l_matannotated, function(x) {
    v <- append(v, as.character(x[,12]))
    return(v)
  })

  col12 <- lapply(col12, function(x) append(x, rep(NA,96 - length(x))))


  new_df <- lapply(seq_along(col12), function(x) {
    mat <- matrix(col12[[x]], nrow = 8, ncol = 12, byrow = TRUE)
    mat[1,8] <- paste("C-",x,sep = "")
    df <- as.data.frame(mat)
    colnames(df) <- 1:12
    rownames(df) <- LETTERS[1:8]
    df %<>% mutate("12" = "LADDER")
    df
  })

  l_matannotated <- lapply(l_matannotated, function(x) {
    x %<>% mutate("12" = "LADDER")
    return(x)
  })

  # Combine with the control names extracted earlier
  lctrls <- list(mol96)
  l_annotated <- c(l_matannotated, lctrls)

  # Append the new "bis" matrices to the list of annotated plates
  l_annotated <- c(l_annotated, new_df)

  return(l_annotated)
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
#' The function reads the specified Excel file, ensuring all required barcode
#' columns are present before proceeding. Barcodes composed from those columns
#' must be non-empty; otherwise the function aborts with a descriptive error.
#' Duplicate barcodes are reported via warnings so users can inspect potential
#' data issues. After validation, the function assembles annotated plate layouts
#' that match the structure of a standard PCR plate.
#' @importFrom readxl read_excel
processFishData <- function(file_path=NULL) {
  file_path <- validate_excel_path(file_path)
  # Import data file and process
  data <- read_excel_safe(file_path)
  validate_required_columns(data, barcode_columns)

  # select features
  data %<>% select(GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP)

  # unite features
  barcode <- data %>% unite(col = BC,
                            GeneName.y,
                            BC1ID,
                            BC1PN,
                            BC1WP,
                            BC2ID,
                            BC2PN,
                            BC2WP,
                            sep = "_")

  validate_barcode_values(barcode$BC, "processFishData() input")

  # select unique rows
  u_barcode <- unique(barcode)

  # number of plates
  nplate <- ceiling(nrow(u_barcode) / 96)

  # replace empty well by NA (last plate)
  u_barcode <- u_barcode$BC
  u_barcode <- append(u_barcode, rep(NA, 96 * nplate - length(u_barcode)))

  # create start and stop sequence to fill plate one by one
  start <- seq(from = 1,
               to = length(u_barcode),
               by = 48)
  stop <- seq(from = 48,
              to = length(u_barcode),
              by = 48)

  # create matrix and fill the matrix with the barcode vector per row
  # ADD C-FLAP and C-
  l_mat <- lapply(seq_along(stop), function(x) {
    mat <-
      matrix(u_barcode[start[x]:stop[x]],
             nrow = 5,
             ncol = 10,
             byrow = TRUE)
    mat[5, 9] <- "C-FLAP"
    mat[5, 10] <- "C-"
    return(mat)
  })

  # annotate row and column of the matrix
  l_matannotated <- lapply(l_mat, function(x) {
    df <- as.data.frame(x)
    colnames(df) <- 2:11
    rownames(df) <- LETTERS[2:6]
    return(df)
  })

  l_matannotated <- lapply(l_matannotated, function(x) {
    x %<>% mutate("1" = NA)
    x %<>% mutate("12" = NA)
    x %<>% select("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
    rownames(x) <- LETTERS[2:6]
    x %<>% rbind(rep(NA, 12))
    x %<>% rbind(rep(NA, 12))
    x %<>% rbind(rep(NA, 12))
    rownames(x)[6] <- "A"
    rownames(x)[7] <- "G"
    rownames(x)[8] <- "H"
    x <- x[order(row.names(x)),]
    return(x)
  })

  #ADD KIF1C & DYNC1H1
  i = 2
  j = 3
  for (z in seq_along(l_matannotated)) {
    l_matannotated[[z]] %<>% mutate_all(as.character)
    rownames(l_matannotated[[z]]) <- LETTERS[1:8]
    l_matannotated[[z]][7, i] <- "KIF1C"
    l_matannotated[[z]][7, j] <- "DYNC1H1"
    i = i + 1
    j = j + 1
    if (i == 11) {
      i = 2
      j = 3
    }
  }

  # get index for the second plate
  idx <- seq(from = 2,
             to = last(seq_along(l_matannotated)),
             by = 2)

  #get name on position 96 and replace value by "C-1"
  v <- vector()
  mol96 <- sapply(seq_along(idx), function(x) {
    v <- append(v, as.character(l_matannotated[[idx[x]]][6, 9]))
    return(v)
  })

  new_l_matannotated <- lapply(seq_along(idx), function(x) {
    l_matannotated[[idx[x]]] %<>% mutate_all(as.character)
    l_matannotated[[idx[x]]][6, 9] <- paste("C-", x, sep = "")
    colnames(l_matannotated[[idx[x]]]) <- 1:12
    rownames(l_matannotated[[idx[x]]]) <- LETTERS[1:8]
    return(l_matannotated[[idx[x]]])
  })

  l_matannotated <- lapply(seq_along(l_matannotated), function(x) {
    if (x  %in% idx)
      l_matannotated[[x]] <- new_l_matannotated[[x / 2]]
    else
      return(l_matannotated[[x]])
  })

  lctrls <- list(mol96)
  l_annotated <- c(l_matannotated, lctrls)
  return(l_annotated)
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
#' The function reads the specified Excel file and validates that all required
#' barcode columns are available before primer fields are stripped. Empty or
#' `NA` barcode values trigger an error, while duplicate barcodes emit a warning
#' to aid troubleshooting. After validation and optional primer removal, the
#' function outputs annotated plate layouts aligned with standard PCR plates.
#' @importFrom readxl read_excel
processFishDataWithoutPrimers <- function(file_path=NULL) {
  file_path <- validate_excel_path(file_path)
  #data <- read_excel(path = args[5])
  data <- read_excel_safe(file_path)
  validate_required_columns(data, barcode_columns)

  # select features
  data %<>% select(GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP)

  # unite features
  barcode <- data %>% unite(col = BC, GeneName.y, BC1ID, BC1PN,
                            BC1WP, BC2ID, BC2PN, BC2WP, sep = "_")

  validate_barcode_values(barcode$BC, "processFishDataWithoutPrimers() input")

  # select unique rows
  u_barcode <- unique(barcode)

  # number of plates
  nplate <- ceiling(nrow(u_barcode)/96)

  # remove primers
  tmp <- lapply(seq_along(u_barcode$BC), function(i) stringi::stri_split_fixed(u_barcode$BC, "_")[[i]][[1]])
  tmp <- do.call("rbind", tmp)
  u_barcode <- as.character(tmp)

  # replace empty well by NA (last plate)
  u_barcode <- append(u_barcode,rep(NA,96*nplate - length(u_barcode)))

  # create start and stop sequence to fill plate one by one
  start <- seq(from = 1, to = length(u_barcode), by = 48)
  stop <- seq(from = 48, to = length(u_barcode), by = 48)

  # create matrix and fill the matrix with the barcode vector per row
  # ADD C-FLAP and C-
  l_mat <- lapply(seq_along(stop), function(x) {
    mat <- matrix(u_barcode[start[x]:stop[x]], nrow = 5, ncol = 10, byrow = TRUE)
    mat[5,9] <- "C-FLAP"; mat[5,10] <- "C-"
    return(mat)
  })

  # annotate row and column of the matrix
  l_matannotated <- lapply(l_mat, function(x) {
    df <- as.data.frame(x)
    colnames(df) <- 2:11
    rownames(df) <- LETTERS[2:6]
    return(df)
  })

  l_matannotated <- lapply(l_matannotated, function(x) {
    x %<>% mutate("1" = NA)
    x %<>% mutate("12" = NA)
    x %<>% select("1","2","3","4","5","6","7","8","9","10","11","12")
    rownames(x) <- LETTERS[2:6]
    x %<>% rbind(rep(NA,12))
    x %<>% rbind(rep(NA,12))
    x %<>% rbind(rep(NA,12))
    rownames(x)[6] <- "A"
    rownames(x)[7] <- "G"
    rownames(x)[8] <- "H"
    x <- x[order(row.names(x)), ]
    return(x)
  })

  #ADD KIF1C & DYNC1H1
  i = 2; j = 3
  for (z in seq_along(l_matannotated)) {
    l_matannotated[[z]] %<>% mutate_all(as.character)
    rownames(l_matannotated[[z]]) <- LETTERS[1:8]
    l_matannotated[[z]][7,i] <- "KIF1C"
    l_matannotated[[z]][7,j] <- "DYNC1H1"
    i = i + 1; j = j + 1
    if (i == 11) {i = 2; j = 3}
  }

  # get index for the second plate
  idx <- seq(from = 2, to =last(seq_along(l_matannotated)), by = 2)

  #get name on position 96 and replace value by "C-1"
  v <- vector()
  mol96 <- sapply(seq_along(idx), function(x){
    v <- append(v, as.character(l_matannotated[[idx[x]]][6,9]))
    return(v)
  })

  new_l_matannotated <- lapply(seq_along(idx), function(x) {
    l_matannotated[[idx[x]]] %<>% mutate_all(as.character)
    l_matannotated[[idx[x]]][6,9] <- paste("C-",x,sep = "")
    colnames(l_matannotated[[idx[x]]]) <- 1:12
    rownames(l_matannotated[[idx[x]]]) <- LETTERS[1:8]
    return(l_matannotated[[idx[x]]])
  })

  l_matannotated <- lapply(seq_along(l_matannotated), function(x) {
    if(x  %in% idx)
      l_matannotated[[x]] <- new_l_matannotated[[x/2]]
    else return(l_matannotated[[x]])
  })

  lctrls <- list(mol96)
  l_annotated <- c(l_matannotated, lctrls)
  return(l_annotated)
}

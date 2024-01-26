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
#' The function reads the user-provided Excel file and selects necessary columns related to the barcode.
#' It then creates a unified barcode by concatenating these columns with an underscore separator.
#' After deduplicating barcodes, it calculates the number of plates needed and arranges barcodes into
#' an 8x12 matrix for each plate. The final list of matrices can be used to guide the setup of PCR plates
#' in a laboratory setting.
#' @importFrom readxl read_excel
getPCR <- function(file_path=NULL) {

  if (is.null(file_path)) {
    # Chemin par défaut vers le fichier inclus dans le package
    file_path <- system.file("extdata", "ResultWithComplementSeq_final_file.xlsx", package = "smFishPlateDesigner")
  }

  if (!file.exists(file_path)) {
    stop("The file doesn't exists: ", file_path)
  }

  # Read the Excel file at the given file path
  data <- readxl::read_excel(path = file_path)

  # Selecting relevant columns for barcode creation
  selected_data <- data %>%
    dplyr::select(GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP)

  # Creating a unified barcode by concatenating selected columns with an underscore
  barcode <- selected_data %>%
    tidyr::unite(col = BC, c(GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP), sep = "_")

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
processDosageTIV <- function(file_path = NULL) {
  if (is.null(file_path)) {
    # Chemin par défaut vers le fichier inclus dans le package
    file_path <- system.file("extdata", "ResultWithComplementSeq_final_file.xlsx", package = "smFishPlateDesigner")
  }

  if (!file.exists(file_path)) {
    stop("The file doesn't exists: ", file_path)
  }
  data <- readxl::read_excel(path = file_path)

  # Select features and create a barcode
  data <- dplyr::select(data, GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP)
  barcode <- tidyr::unite(data, col = BC, GeneName.y, BC1ID, BC1PN,
                          BC1WP, BC2ID, BC2PN, BC2WP, sep = "_")

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
#' The function reads the specified Excel file, extracting necessary information for barcode creation
#' and plate layout organization. It then processes this information to output a list of plate layouts,
#' each corresponding to a unique set of samples and annotations.
#' @importFrom readxl read_excel
processFishData <- function(file_path=NULL) {
  if (is.null(file_path)) {
    # Chemin par défaut vers le fichier inclus dans le package
    file_path <- system.file("extdata", "ResultWithComplementSeq_final_file.xlsx", package = "smFishPlateDesigner")
  }

  if (!file.exists(file_path)) {
    stop("The file doesn't exists: ", file_path)
  }
  # Import data file and process
  data <- readxl::read_excel(path = file_path)

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
#' The function reads the specified Excel file, extracting necessary information for barcode creation
#' and plate layout organization, excluding primer information. It processes this information to output
#' a list of plate layouts, each corresponding to a unique set of samples and annotations.
#' @importFrom readxl read_excel
processFishDataWithoutPrimers <- function(file_path=NULL) {
  if (is.null(file_path)) {
    # Chemin par défaut vers le fichier inclus dans le package
    file_path <- system.file("extdata", "ResultWithComplementSeq_final_file.xlsx", package = "smFishPlateDesigner")
  }

  if (!file.exists(file_path)) {
    stop("The file doesn't exists: ", file_path)
  }
  #data <- read_excel(path = args[5])
  data <- readxl::read_excel(path = file_path)

  # select features
  data %<>% select(GeneName.y, BC1ID, BC1PN, BC1WP, BC2ID, BC2PN, BC2WP)

  # unite features
  barcode <- data %>% unite(col = BC, GeneName.y, BC1ID, BC1PN,
                            BC1WP, BC2ID, BC2PN, BC2WP, sep = "_")

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

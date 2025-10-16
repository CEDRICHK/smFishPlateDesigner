#' Probes dataset used in examples
#'
#' Structure and column metadata for the Excel file distributed under
#' `inst/extdata/Probes_final_Oligopool11_with_sets.xlsx`. The file contains
#' 34 columns describing HT-smFISH probe design results.
#'
#' @details
#' Columns include identifiers (`ENST`, `ENSG.x`, `GeneName.x`), probe metrics
#' (`dGOpt`, `ProbeSize`, `Seq`, etc.), quality filters, and barcode fields
#' (`BC1*`, `BC2*`, `WithComplementSeq`). Users should load the Excel file with
#' `readxl::read_excel()` when they need the actual data content.
#' @return Invisibly returns a character vector of column names.
#' @export
probes_final_columns <- function() {
  c(
    "ENST", "ENSG.x", "GeneName.x", "SET", "ENSG.y", "GeneName.y",
    "dGOpt", "theStartPos", "theEndPos", "ProbeSize", "Seq", "dGScore",
    "dG37", "GCpc", "GCFilter", "aCompFilter", "aStackFilter",
    "cCompFilter", "cStackFilter", "cSpecStackFilter", "NbOfPNAS",
    "PNASFilter", "RSESeqFilter", "InsideUTR", "BC1ID", "BC1PN",
    "BC1WP", "BC1", "BC2ID", "BC2PN", "BC2WP", "BC2",
    "BC1YHybXBC2", "WithComplementSeq"
  )
}

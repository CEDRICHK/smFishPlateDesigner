#' HT-smFish Result Dataset: ResultWithComplementSeq_final_file
#'
#' This dataset contains detailed results for  High-Throughput single-molecule fluorescence in situ hybridization (HT-smFish) analyses, including plate layouts and complementary sequences. The dataset encompasses a range of metrics and identifiers relevant to the HT-smFish experiments.
#'
#' @format A `data.frame` with X rows and 33 columns.
#' @details
#'   - `ENST`: Transcript identifier.
#'   - `ENSG.x`: Gene identifier (x-version).
#'   - `GeneName.x`: Gene name (x-version).
#'   - `SET`: Experimental set identifier.
#'   - `ENSG.y`: Gene identifier (y-version).
#'   - `GeneName.y`: Gene name (y-version).
#'   - `dGOpt`: Optimal delta G value.
#'   - `theStartPos`: Start position of the sequence.
#'   - `theEndPos`: End position of the sequence.
#'   - `ProbeSize`: Size of the probe.
#'   - `Seq`: Sequence of nucleotides.
#'   - `dGScore`: Delta G score.
#'   - `dG37`: Delta G value at 37 degrees Celsius.
#'   - `GCpc`: GC content percentage.
#'   - `GCFilter`, `aCompFilter`, `aStackFilter`, `cCompFilter`, `cStackFilter`, `cSpecStackFilter`: Various sequence filters.
#'   - `NbOfPNAS`: Number of PNAS.
#'   - `PNASFilter`: PNAS filter status.
#'   - `RSESeqFilter`: RSE sequence filter status.
#'   - `InsideUTR`: Indicator if inside UTR.
#'   - `BC1ID`, `BC1PN`, `BC1WP`, `BC1`: Barcode 1 details.
#'   - `BC2ID`, `BC2PN`, `BC2WP`, `BC2`: Barcode 2 details.
#'   - `BC1YHybXBC2`: Hybrid barcode details.
#'   - `WithComplementSeq`: Complementary sequence information.
#'
#' @source Describe the source of the data if applicable.
"ResultWithComplementSeq_final_file"

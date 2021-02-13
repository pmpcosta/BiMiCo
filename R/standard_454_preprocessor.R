## TODO: define optional args for read-filtering
##        input filename not only "fastq"
##        output dir not only subdir of input dir, overwrite permission

#' standard preprocessor for 454 reads in BiMiCo pipeline
#'
#' Convenience wrapper for DADA2 read processing functions. Given a directory of raw 454 fastq files, performs platform-specific trimming and QC, writes reads to output directory.
#' @param readfiles (Required) Path to single-end 454 demultiplexed fastq files
#' @param outdir (Required) Path to write quality-filtered fastq files to
#' @param trim_read_length (Optional) max. length at which to truncate reads. Default = 0 (no truncating)
#' @param trim5end (Optional) trim n reads from the start of reads (unidentified adapters, barcodes). Default=0 (no trimming)
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio). Default=F
#' @keywords read processing dada2
#' @export
#' @examples
#' prep_454()

prep_454 <- function(readfiles, outdir, mtthread=F, trim_read_length=0, trim5end=0){

  # store fastq file names
  fqs <- sort(
    list.files(
      readfiles,
      full.names = TRUE
      )
    )

  # extract sample names
  samps <- sapply(
    strsplit(
    basename(fqs),
    ".fastq"),
    `[`,
    1
    )

  # create dir for filtered, trimmed fastq files
  filt_fqs <- file.path(
    outdir,
    paste0(samps, "_filt.fastq.gz")
    )


  # default filter parameters to work with DADA2 (454)
  fout <- dada2::filterAndTrim(
    fqs,
    filt_fqs,
                        maxN=0,
                        maxEE=c(2),
                        truncQ=2,
                        truncLen = trim_read_length,
                        trimLeft = trim5end,
                        rm.phix=TRUE,
                        compress=TRUE,
                        multithread=mtthread
    )
# Output filter summary

  closeAllConnections()

  return(fout)

  }

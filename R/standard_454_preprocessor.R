## TODO: define optional args for read-filtering
##        input filename not only "fastq"
##        output dir not only subdir of input dir, overwrite permission

#' standard preprocessor for 454 reads in BiMiCo pipeline
#'
#' Convenience wrapper for DADA2 read processing functions. Given a directory of raw 454 fastq files, performs platform-specific trimming and QC, outputs reads to a sub-folder of the specified name
#' @param readfiles (Required) Path to single-end 454 demultiplexed fastq files
#' @param outdir (Required) String specifying output directory for processed files; will be created as subdirectory of input dir.
#' @param trim_read_length (Optional) max. length at which to truncate reads. Default = 0 (no truncating)
#' @param mtthread (Required) Boolean, enables multithreading (not recommended in Rstudio)
#' @keywords read processing dada2
#' @export
#' @examples
#' std_454_preprocess()

std_454_preprocess <- function(readfiles, outdir, mtthread, trim_read_length=0){

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
    readfiles,
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
                        rm.phix=TRUE,
                        compress=TRUE,
                        multithread=mtthread
    )
# Output filter summary

  closeAllConnections()

  return(fout)

  }

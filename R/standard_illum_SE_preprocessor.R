## TODO: define optional args for read-filtering
##        input filename not only "fastq"
##        output dir not only subdir of input dir, overwrite permission

#' standard preprocessor for Illumina single-end reads in BiMiCo pipeline
#'
#' UNDER CONSTRUCTION
#' @param readfiles (Required) Path to single-end Illumina demultiplexed fastq files
#' @param outdir (Required) String specifying output directory for processed files; will be created as subdirectory of input dir.
#' @param trim_read_length (Optional) max. length at which to truncate reads. Default = 0 (no truncating)
#' @keywords read processing dada2
#' @export
#' @examples
#' std_illum_SE_preprocess()

std_illum_SE_preprocess <- function(readfiles, outdir, trim_read_length=0){

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
    p1,
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
                        multithread=FALSE
    )
# Output filter summary

  closeAllConnections()

  return(fout)

  }

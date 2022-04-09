#' @export
#' @name loadSamplesFile
#' @title Loading sample file
#' @author Haniel Cedraz
#' @details December 2021
#' @usage
#' loadSamplesFile(file, reads_folder, column = "SAMPLE_ID", libraryType = "pairEnd")
#' @description
#' Function to load the sample file
#'
#' @param file
#' \code{Character.} The filename of the sample file. Default samples.txt.
#' @param reads_folder
#' \code{Character.} Directory where the raw sequence data is stored. Default 00-Fastq.
#' @param column
#' \code{Character.} Column name from the sample sheet to use as read folder names. Default SAMPLE_ID
#' @param libraryType
#' \code{Character.} The library type to use. Available: 'pairEnd' or 'singleEnd'. Default pairEnd
#' @importFrom tools file_ext
#' @importFrom dplyr mutate %>%
#' @importFrom glue glue
#' @importFrom readr read_table
#' @importFrom purrr modify_if
#' @export



## loadSampleFile
loadSamplesFile <- function(file, reads_folder, column = "SAMPLE_ID", libraryType = "pairEnd"){
  ## debug
  #file = opt$samplesFile; reads_folder = opt$Raw_Folder; column = opt$samplesColumn
  ##

  aceptedLibraryType <- c("pairEnd", "singleEnd")
  if (!libraryType %in% aceptedLibraryType) {
    stop(glue("Library type ({libraryType}) not found, please provide one of 'pairEnd' or 'singleEnd'"))
  }

  if (!file.exists(file)) {
    stop(glue("Sample file {file} does not exist"))
  }
  ### column SAMPLE_ID should be the sample name
  ### rows can be commented out with #
  if (libraryType == "singleEnd") {
    targets <- read_table(file, col_types = list("c", "c")) %>%
      modify_if(~is.double(.), ~as.character(.)) %>%
      as.data.frame()

  } else if (libraryType == "pairEnd") {
    targets <- read_table(file, col_types = list("c", "c", "c")) %>%
      modify_if(~is.double(.), ~as.character(.)) %>%
      as.data.frame()
  }


  if (libraryType == "pairEnd") {
    if (!all(c(column, "Read_1", "Read_2") %in% colnames(targets))) {
      stop(glue("Expecting the three columns SAMPLE_ID, Read_1 and Read_2 in samples file (tab-delimited)"))
      stop()
    }
  }
  #for (i in seq.int(nrow(targets$column))) {
    #if (targets[i, column]) {
      ext <- unique(file_ext(dir(file.path(reads_folder), pattern = "fastq|gz")))
      if (length(ext) == 0) {
        write(paste("Cannot locate fastq or sff file in folder", reads_folder, "\n"), stderr())
        stop()
      }
      # targets$type[i] <- paste(ext,sep="/")
    #}
    # else {
    #   ext <- file_ext(grep("gz", dir(file.path(reads_folder)), value = TRUE))
    #   if (length(ext) == 0) {
    #     write(paste(targets[i,column],"is not a gz file\n"), stderr())
    #     stop()
    #   }
    #
    # }
  #}
  write(glue("{file} contains {nrow(targets)} samples to process"), stdout())
  return(targets)
}

# library(dplyr)
# library(glue)
# library(tools)
# library(purrr)
# library(readr)
# samplesFile <- "/Users/haniel/OneDrive/posDoc/miRbaqcom/samplesReal.txt"
# cleanedFolder <- "/Users/haniel/OneDrive/posDoc/miRbaqcom/01-CleanedReadsReal/"
# samplesColumn <- "SAMPLE_ID"
# slibraryType <- "singleEnd"
# samples <- loadSamplesFile(
#   file = samplesFile,
#   reads_folder = cleanedFolder,
#   column = samplesColumn,
#   libraryType = slibraryType
# )

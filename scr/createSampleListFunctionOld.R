#' @export
#' @name createSampleListOld
#' @title Creating sample list to run mcapply parallel
#' @author Haniel Cedraz
#' @details December 2021
#' @usage
#' createSampleListOld(samples, reads_folder, column = "SAMPLE_ID", fileType = NULL, libraryType = "pairEnd", samplesFromSTAR = FALSE, step = NULL)
#' @description
#' Function to create sample list to run mcapply parallel
#'
#' @param samples
#' \code{Character.} The filename of the sample file. Default samples.txt.
#' @param reads_folder
#' \code{Character.} Directory where the raw sequence data is stored. Default 00-Fastq.
#' @param column
#' \code{Character.} Column name from the sample sheet to use as read folder names. Default SAMPLE_ID
#' @param fileType
#' \code{Character.} The file type to use. Available: 'fastq.gz', 'bam' or 'sam'. Default NULL
#' @param libraryType
#' \code{Character.} The library type to use. Available: 'pairEnd' or 'singleEnd'. Default pairEnd
#' @param samplesFromSTAR
#' \code{logical.} Whether want to count reads from STAR mapped files. Default FALSE
#' @param step
#' \code{Character} Which step to run, Quality Controls or Mapping. Default NULL
#' @importFrom glue glue
#' @importFrom utils read.table
#' @export




createSampleListOld <- function(samples, reads_folder, column = "SAMPLE_ID", fileType = NULL, libraryType = "pairEnd", samplesFromSTAR = FALSE, step = NULL) {

  sampleList <- list()
  samples <- as.data.frame(samples)
  aceptedFileTypes <- c("fastq.gz", "bam", "sam")


if (!is.null(fileType)) {
  if (!fileType %in% aceptedFileTypes) {
    stop(glue("File type ({fileType}) not found, please provide one of 'fastq', 'bam' or 'sam'"))
  }
}


  aceptedLibraryType <- c("pairEnd", "singleEnd")
  if (!libraryType %in% aceptedLibraryType) {
    stop(glue("Library type ({libraryType}) not found, please provide one of 'pairEnd' or 'singleEnd'"))
  }


  if (!file.exists(reads_folder)) {
    stop(glue("reads_folder {reads_folder} does not exist\n"))

  }
  ### column SAMPLE_ID should be the sample name
  ### rows can be commented out with #
  if (libraryType == "pairEnd") {
    if (!all(c(column, "Read_1", "Read_2") %in% colnames(samples))) {
      stop(glue("Expecting the three columns {column}, Read_1 and Read_2 in samples file (tab-delimited)\n"))

    }
  }

if (!is.null(fileType)) {
  if (fileType == "fastq.gz") {
    if (libraryType == "pairEnd") {
      sampleList <- list()
      for (i in 1:nrow(samples)) {
        reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
        #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$")
        if (step == "QualityControl") {
          map <- lapply(c("_R1","_R2"), grep, x = reads, value = TRUE)
          names(map) <- c("R1","R2")
          map$sampleName <-  samples[i,column]
          map$R1 <- samples[i,2]
          map$R2 <- samples[i,3]
          sampleList[[paste(map$sampleName)]] <- map
          #sampleList[[paste(map$sampleName)]]
        } else if (step == "Mapping") {
          map <- lapply(c("_PE1", "_PE2", "_SE1", "_SE2"), grep, x = reads, value = TRUE)
          names(map) <- c("PE1", "PE2", "SE1", "SE2")
          map$sampleName <-  samples[i,column]
          map$PE1 <- map$PE1[i]
          map$PE2 <- map$PE2[i]
          map$SE1 <- map$SE1[i]
          map$SE2 <- map$SE2[i]
          sampleList[[paste(map$sampleName)]] <- map
          sampleList[[paste(map$sampleName, sep = "_")]]
        }

      }
      write(paste("Setting up",length(sampleList),"jobs"), stdout())
      return(sampleList)

    } else if (libraryType == "singleEnd") {

      if (step == "bowtie") {

        sampleList <- list()
        for (i in 1:nrow(samples)) {
          #reads_folder <- "01-CleanedReads/"
          reads <- dir(path = file.path(reads_folder), pattern = "fastq$")
          #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$")
          map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
          names(map) <- c("SE")
          map$sampleName <-  samples[i,column]
          map$SE <- map$SE[i]
          #map$SE <- samples[i,2]
          sampleList[[paste(map$sampleName)]] <- map
          #sampleList[[paste(map$sampleName)]]

        }
        write(paste("Setting up",length(sampleList),"jobs"), stdout())
        return(sampleList)


      } else {
        sampleList <- list()
        for (i in 1:nrow(samples)) {
          #reads_folder <- "01-CleanedReads/"
          reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
          #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$")
          if (any(step %in% c("QualityControl", "Mapping"))) {
            map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
            names(map) <- c("SE")
            map$sampleName <-  samples[i,column]
            map$SE <- map$SE[i]
            #map$SE <- samples[i,2]
            sampleList[[paste(map$sampleName)]] <- map
            #sampleList[[paste(map$sampleName)]]
          }

        }
        write(paste("Setting up",length(sampleList),"jobs"), stdout())
        return(sampleList)
      }

    }

  } else if (fileType == "bam") {
    sampleList <- list()
    for (i in 1:nrow(samples)) {
      reads <- dir(path = file.path(reads_folder), pattern = "bam$", full.names = TRUE)
      map <- lapply(c("_sam_sorted_pos.bam"), grep, x = reads, value = TRUE)
      names(map) <- c("bam_sorted_pos")
      map$sampleName <-  samples[i,column]
      map$bam_sorted_pos <- map$bam_sorted_pos[i]

      sampleList[[paste(map$sampleName)]] <- map
      sampleList[[paste(map$sampleName, sep = "_")]]
    }
    write(paste("Setting up",length(sampleList),"jobs"), stdout())
    return(sampleList)

  } else if (fileType == "sam") {
    sampleList <- list()
    for (i in 1:nrow(samples)) {
      reads <- dir(path = file.path(reads_folder), pattern = "sam$", full.names = TRUE)
      map <- lapply(c("_unsorted_sample.sam"), grep, x = reads, value = TRUE)
      names(map) <- c("unsorted_sample")
      map$sampleName <-  samples[i,column]
      map$unsorted_sample <- map$unsorted_sample[i]

      sampleList[[paste(map$sampleName)]] <- map
      sampleList[[paste(map$sampleName, sep = "_")]]
    }
    write(paste("Setting up",length(sampleList),"jobs"), stdout())
    return(sampleList)
  }
}

  if (is.null(fileType) & samplesFromSTAR == TRUE) {
    # if (is.null(starFolder)) {
    #   stop(glue("{starFolder} cannot be null, please provide a valid name for the STAR results folder"))
    # }
    #reads_folder <- starFolder
    sampleList <- list()
      for (i in 1:nrow(samples)) {
        reads <- dir(path = file.path(reads_folder), pattern = ".bam$", full.names = TRUE)

        map <- lapply(c("_Aligned.sortedByCoord.out.bam"), grep, x = reads, value = TRUE)
        names(map) <- c("Aligned.sortedByCoord.out")
        map$sampleName <-  samples[i,column]
        map$Aligned.sortedByCoord.out <- map$Aligned.sortedByCoord.out[i]

        sampleList[[paste(map$sampleName)]] <- map
        sampleList[[paste(map$sampleName, sep = "_")]]
      }
      write(paste("Setting up", length(sampleList), "jobs"),stdout())
      return(sampleList)
  }

  return(sampleList)
}

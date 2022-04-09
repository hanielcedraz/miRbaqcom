#' @export
#' @name createSampleList
#' @title Creating sample list to run mcapply parallel
#' @author Haniel Cedraz
#' @details December 2021
#' @usage
#' createSampleList(samples, reads_folder, column, program, libraryType, fileType, fromSTAR = FALSE)
#' @description
#' Function to create sample list to run mcapply parallel
#'
#' @param samples
#' \code{Character.} The filename of the sample file. Default samples.txt.
#' @param reads_folder
#' \code{Character.} Directory where the raw sequence data is stored. Default 00-Fastq.
#' @param column
#' \code{Character.} Column name from the sample sheet to use as read folder names. Default SAMPLE_ID
#' @param program
#' \code{Character.} The program to be used. Available 'trimmomatic', 'star', 'Hisat2', 'htseq', 'featurecount', 'bowtie', 'samtools' or 'bwa'
#' @param libraryType
#' \code{Character.} The library type to use. Available: 'pairEnd' or 'singleEnd'.
#' @param fileType
#' \code{Character.} The file type to use. Available: 'bam' or 'sam'.
#' @param fromSTAR
#' \code{Logical.} Wheter to use samples from STAR in FeatureCouts
#' @importFrom glue glue
#' @export



#reads_folder <- "/Users/haniel/OneDrive/posDoc/miRbaqcom/01-CleanedReadsReal/"
createSampleList <- function(samples, reads_folder, column, program, libraryType, fileType = NULL, fromSTAR = FALSE) {

  samples <- as.data.frame(samples)
  aceptedFileTypes <- c("bam", "sam", "fastq.gz")

  if (!is.null(fileType)){
    if (!fileType %in% aceptedFileTypes) {
      stop(glue("File type ({fileType}) not found, please provide one of 'bam' or 'sam' "))
    }
  }


  aceptedPrograms <- c('trimmomatic', 'star', 'Hisat2', 'htseq', 'featurecount', 'bowtie', 'samtools', 'bwa')
  if (!program %in% aceptedPrograms) {
    stop(glue("Program ({program}) not found, please provide one of 'trimmomatic', 'star', 'Hisat2', 'htseq', 'featurecount', 'bowtie', 'samtools' or 'bwa'"))
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




  if (program == "trimmomatic") {
    if (libraryType == "pairEnd") {
      #qcList <- function(samples, reads_folder, column){
      mapping_list <- list()
      for (i in 1:nrow(samples)) {
        reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
        #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
        map <- lapply(c("_R1","_R2"), grep, x = reads, value = TRUE)
        names(map) <- c("R1","R2")
        map$sampleName <-  samples[i,column]
        map$R1 <- samples[i,2]
        map$R2 <- samples[i,3]
        mapping_list[[paste(map$sampleName)]] <- map
        #mapping_list[[paste(map$sampleName)]]
      }
      write(paste("Setting up",length(mapping_list),"jobs"), stdout())
      return(mapping_list)
      #}
    } else if (libraryType == "singleEnd") {
      #qcList <- function(samples, reads_folder, column){
      mapping_list <- list()
      # samples <- read.table("samplesSingle.txt", header = TRUE)
      # reads_folder <- "00-FastqSingle/"
      # column <- "SAMPLE_ID"
      for (i in 1:nrow(samples)) {

        reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
        #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
        map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
        names(map) <- c("SE")
        map$sampleName <-  samples[i,column]
        #map$SE <- map$SE[i]
        map$SE <- map$SE[i]
        mapping_list[[paste(map$sampleName)]] <- map
        #mapping_list[[paste(map$sampleName)]]
      }
      write(paste("Setting up",length(mapping_list),"jobs"), stdout())
      return(mapping_list)
      #}
    }
  } else if (program == "star" | program == "Hisat2" | program == "bowtie" | program == "bwa") {
    if (libraryType == "pairEnd") {
      #mappingList <- function(samples, reads_folder, column){
      mapping_list <- list()

      for (i in 1:nrow(samples)) {
        reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
        if (program == "bowtie") {
          reads <- dir(path = file.path(reads_folder), pattern = "fastq$", full.names = TRUE)
        }
        # for (i in seq.int(to=nrow(samples))){
        #     reads <- dir(path=file.path(reads_folder,samples[i,column]),pattern="gz$",full.names=TRUE)
        map <- lapply(c("_PE1", "_PE2", "_SE1", "_SE2"), grep, x = reads, value = TRUE)
        names(map) <- c("PE1", "PE2", "SE1", "SE2")
        map$sampleName <-  samples[i,column]
        map$PE1 <- map$PE1[i]
        map$PE2 <- map$PE2[i]
        map$SE1 <- map$SE1[i]
        map$SE2 <- map$SE2[i]
        mapping_list[[paste(map$sampleName)]] <- map
        mapping_list[[paste(map$sampleName, sep = "_")]]
      }
      write(paste("Setting up", length(mapping_list), "jobs"),stdout())
      return(mapping_list)
      #}

    } else if (libraryType == "singleEnd") {
      #mappingList <- function(samples, reads_folder, column){
      mapping_list <- list()
      for (i in 1:nrow(samples)) {
        #reads_folder <- "01-CleanedReadsSingle/"
        reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
        if (program == "bowtie") {
          reads <- dir(path = file.path(reads_folder), pattern = "fastq$", full.names = TRUE)
        }
        #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
        map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
        names(map) <- c("SE")
        map$sampleName <-  samples[i,column]
        map$SE <- map$SE[i]
        #map$R2 <- samples[i,3]
        mapping_list[[paste(map$sampleName)]] <- map
        #mapping_list[[paste(map$sampleName)]]
      }
      write(paste("Setting up",length(mapping_list),"jobs"), stdout())
      return(mapping_list)
      #}
    }
  } else if (program == "htseq") {
    counting_list <- list()
    if (casefold(fileType, upper = FALSE) == 'bam') {
      for (i in 1:nrow(samples)) {
        files <- dir(path = file.path(reads_folder), recursive = TRUE, pattern = paste0('.bam$'), full.names = TRUE)

        count <- lapply(c("_sam_sorted_pos.bam"), grep, x = files, value = TRUE)
        names(count) <- c("bam_sorted_pos")
        count$sampleName <-  samples[i,column]
        count$bam_sorted_pos <- count$bam_sorted_pos[i]

        counting_list[[paste(count$sampleName)]] <- count
        counting_list[[paste(count$sampleName, sep = "_")]]

      }
    } else if (casefold(fileType, upper = FALSE) == 'sam') {
      for (i in 1:nrow(samples)) {
        files <- dir(path = file.path(reads_folder), recursive = TRUE, pattern = paste0('.sam$'), full.names = TRUE)

        count <- lapply(c("_unsorted_sample.sam"), grep, x = files, value = TRUE)
        names(count) <- c("unsorted_sample")
        count$sampleName <-  samples[i,column]
        count$unsorted_sample <- count$unsorted_sample[i]

        counting_list[[paste(count$sampleName)]] <- count
        counting_list[[paste(count$sampleName, sep = "_")]]

      }
    }
    write(paste("Setting up", length(counting_list), "jobs"),stdout())
    return(counting_list)

  } else if (program == "featurecount") {
    if (fromSTAR) {
      #reads_folder <- opt$inputFolder
      #countingList <- function(samples, reads_folder, column) {
      counting_list <- list()
      for (i in 1:nrow(samples)) {
        files <- dir(path = file.path(reads_folder), recursive = TRUE, pattern = paste0('.bam$'), full.names = TRUE)

        count <- lapply(c("_Aligned.sortedByCoord.out.bam"), grep, x = files, value = TRUE)
        names(count) <- c("Aligned.sortedByCoord.out")
        count$sampleName <-  samples[i,column]
        count$Aligned.sortedByCoord.out <- count$Aligned.sortedByCoord.out[i]

        counting_list[[paste(count$sampleName)]] <- count
        counting_list[[paste(count$sampleName, sep = "_")]]
      }
      write(paste("Setting up", length(counting_list), "jobs"),stdout())
      return(counting_list)
      #}

    } else {
      #countingList <- function(samples, reads_folder, column) {
      counting_list <- list()
      if (casefold(fileType, upper = FALSE) == 'bam') {
        for (i in 1:nrow(samples)) {
          files <- dir(path = file.path(reads_folder), recursive = TRUE, pattern = paste0('.bam$'), full.names = TRUE)

          count <- lapply(c("_sam_sorted_pos.bam"), grep, x = files, value = TRUE)
          names(count) <- c("bam_sorted_pos")
          count$sampleName <-  samples[i,column]
          count$bam_sorted_pos <- count$bam_sorted_pos[i]

          counting_list[[paste(count$sampleName)]] <- count
          counting_list[[paste(count$sampleName, sep = "_")]]

        }
      } else if (casefold(fileType, upper = FALSE) == 'sam') {
        for (i in 1:nrow(samples)) {
          files <- dir(path = file.path(reads_folder), recursive = TRUE, pattern = paste0('.sam$'), full.names = TRUE)

          count <- lapply(c("_unsorted_sample.sam"), grep, x = files, value = TRUE)
          names(count) <- c("unsorted_sample")
          count$sampleName <-  samples[i,column]
          count$unsorted_sample <- count$unsorted_sample[i]

          counting_list[[paste(count$sampleName)]] <- count
          counting_list[[paste(count$sampleName, sep = "_")]]

        }
      }
      write(paste("Setting up", length(counting_list), "jobs"),stdout())
      return(counting_list)
      #}
    }
  } else if (program == "samtools") {
    #samtoolsList <- function(samples, reads_folder, column){
    samtoolsfiles <- list()
    for (i in 1:nrow(samples)) {
      samfiles <- dir(path = file.path(reads_folder), recursive = TRUE, pattern = ".sam$", full.names = TRUE)
      maps <- lapply(c("_unsorted_sample"), grep, x = samfiles, value = TRUE)
      names(maps) <- c("unsorted_sample")
      maps$sampleName <-  samples[i,column]
      maps$unsorted_sample <- maps$unsorted_sample[i]

      samtoolsfiles[[paste(maps$sampleName)]] <- maps
      samtoolsfiles[[paste(maps$sampleName, sep = "_")]]

    }
    write(paste("Setting up", length(samtoolsfiles), "jobs"),stdout())
    return(samtoolsfiles)
    #}

  }


}




# else if (program == "bowtie") {
#     if (libraryType == "pairEnd") {
#         sampleList <- list()
#         for (i in 1:nrow(samples)) {
#             reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
#             #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
#             map <- lapply(c("_R1","_R2"), grep, x = reads, value = TRUE)
#             names(map) <- c("R1","R2")
#             map$sampleName <-  samples[i,column]
#             map$R1 <- samples[i,2]
#             map$R2 <- samples[i,3]
#             sampleList[[paste(map$sampleName)]] <- map
#             #sampleList[[paste(map$sampleName)]]
#         }
#         write(paste("Setting up",length(sampleList),"jobs"), stdout())
#         return(sampleList)
#     } else if (libraryType == "singleEnd") {
#         sampleList <- list()
#         for (i in 1:nrow(samples)) {
#             reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
#             #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
#             map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
#             names(map) <- c("SE")
#             map$sampleName <-  samples[i,column]
#             #map$SE <- map$SE[i]
#             map$SE <- samples[i,2]
#             sampleList[[paste(map$sampleName)]] <- map
#             #sampleList[[paste(map$sampleName)]]
#         }
#         write(paste("Setting up",length(sampleList),"jobs"), stdout())
#         return(sampleList)
#     }
# }   # else if (program == "Hisat2") {
#     if (libraryType == "pairEnd") {
#         mappingList <- function(samples, reads_folder, column){
#             mapping_list <- list()
#             for (i in 1:nrow(samples)) {
#                 reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
#                 # for (i in seq.int(to=nrow(samples))){
#                 #     reads <- dir(path=file.path(reads_folder,samples[i,column]),pattern="gz$",full.names=TRUE)
#                 map <- lapply(c("_PE1", "_PE2", "_SE1", "_SE2"), grep, x = reads, value = TRUE)
#                 names(map) <- c("PE1", "PE2", "SE1", "SE2")
#                 map$sampleName <-  samples[i,column]
#                 map$PE1 <- map$PE1[i]
#                 map$PE2 <- map$PE2[i]
#                 map$SE1 <- map$SE1[i]
#                 map$SE2 <- map$SE2[i]
#                 mapping_list[[paste(map$sampleName)]] <- map
#                 mapping_list[[paste(map$sampleName, sep = "_")]]
#             }
#             write(paste("Setting up", length(mapping_list), "jobs"),stdout())
#             return(mapping_list)
#         }
#
#     } else if (libraryType == "singleEnd") {
#         mappingList <- function(samples, reads_folder, column){
#             mapping_list <- list()
#             for (i in 1:nrow(samples)) {
#                 reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
#                 #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
#                 map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
#                 names(map) <- c("SE")
#                 map$sampleName <-  samples[i,column]
#                 map$SE <- map$SE[i]
#                 #map$R2 <- samples[i,3]
#                 mapping_list[[paste(map$sampleName)]] <- map
#                 #mapping_list[[paste(map$sampleName)]]
#             }
#             write(paste("Setting up",length(mapping_list),"jobs"), stdout())
#             return(mapping_list)
#         }
#     }
# }




#
#
#
#
#
#
#
#
#
# if (fileType == "fastq.gz") {
#     if (libraryType == "pairEnd") {
#         sampleList <- list()
#         for (i in 1:nrow(samples)) {
#             reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
#             #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
#             map <- lapply(c("_R1","_R2"), grep, x = reads, value = TRUE)
#             names(map) <- c("R1","R2")
#             map$sampleName <-  samples[i,column]
#             map$R1 <- samples[i,2]
#             map$R2 <- samples[i,3]
#             sampleList[[paste(map$sampleName)]] <- map
#             #sampleList[[paste(map$sampleName)]]
#         }
#         write(paste("Setting up",length(sampleList),"jobs"), stdout())
#         return(sampleList)
#     } else if (libraryType == "singleEnd") {
#         sampleList <- list()
#         for (i in 1:nrow(samples)) {
#             reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
#             #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
#             map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
#             names(map) <- c("SE")
#             map$sampleName <-  samples[i,column]
#             #map$SE <- map$SE[i]
#             map$SE <- samples[i,2]
#             sampleList[[paste(map$sampleName)]] <- map
#             #sampleList[[paste(map$sampleName)]]
#         }
#         write(paste("Setting up",length(sampleList),"jobs"), stdout())
#         return(sampleList)
#     }
#
# } else if (fileType == "bam") {
#     sampleList <- list()
#     for (i in 1:nrow(samples)) {
#         reads <- dir(path = file.path(reads_folder), pattern = "bam$", full.names = TRUE)
#         map <- lapply(c(".bam"), grep, x = reads, value = TRUE)
#         names(map) <- c("bam")
#         map$sampleName <-  samples[i,column]
#         map$bam_sorted_pos <- map$bam[i]
#
#         sampleList[[paste(map$sampleName)]] <- map
#         sampleList[[paste(map$sampleName, sep = "_")]]
#     }
#     write(paste("Setting up",length(sampleList),"jobs"), stdout())
#     return(sampleList)
#
# } else if (fileType == "sam") {
#     sampleList <- list()
#     for (i in 1:nrow(samples)) {
#         reads <- dir(path = file.path(reads_folder), pattern = "sam$", full.names = TRUE)
#         map <- lapply(c(".sam"), grep, x = reads, value = TRUE)
#         names(map) <- c("sam")
#         map$sampleName <-  samples[i,column]
#         map$bam_sorted_pos <- map$bam[i]
#
#         sampleList[[paste(map$sampleName)]] <- map
#         sampleList[[paste(map$sampleName, sep = "_")]]
#     }
#     write(paste("Setting up",length(sampleList),"jobs"), stdout())
#     return(sampleList)
# }





#}

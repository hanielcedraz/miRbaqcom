#!/usr/bin/env Rscript


########################################
### LOADING PACKAGES
########################################
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("baqcomPackage"))

########################################
### SETING PARAMETERS
########################################
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")


option_list <- list(
    # make_option(
    #     opt_str = c("-i", "--install"),
    #     type = "logical",
    #     default = TRUE,
    #     help = "Wether install UMI-Tools or not",
    #     dest = "installUMItools"
    # ),
    make_option(
        opt_str = c("-f", "--file"), 
        type = "character", 
        default = "samples.txt",
        help = "The filename of the sample file [default %default]",
        dest = "samplesFile"
    ),
    make_option(
        opt_str = c("-i", "--inputDirectory"), 
        type = "character", 
        default = "00-Fastq",
        help = "Directory where the raw sequence data is stored [default %default]",
        dest = "rawFolder"
    ),
    make_option(
        opt_str = c("-c", "--column"), 
        type = "character", 
        default = "SAMPLE_ID",
        help = "Column name from the sample sheet to use as read folder names [default %default]",
        dest = "samplesColumn"
    ),
    make_option(
        opt_str = c("-O", "--umiOption"),
        type = "character",
        default = "extract",
        help = glue(
            "umi_tools command to use. Options: 'extract', 'dedup' or 'both'. [default %default]."
        ),
        dest = "umiCommand"
    ),
    make_option(
        opt_str = c("-e", "--extractMethod"),
        type = "character",
        default = "string",
        help = glue("If %umiOption extract, how to extract the umi +/- cell barcodes, Choose from 'string' or 'regex' [default %default]"),
        dest = "extractMethod"
    ),
    make_option(
        opt_str = c("-d", "--discardRegex"),
        action = 'store_true',
        type = "logical",
        default = FALSE,
        help = "If %extractMethod regex, wheter to discard bases before and after barcode/umi  [default %default]",
        dest = "discard"
    ),
    make_option(
        opt_str = c("-a", "--bcPattern"), 
        type  = 'character', 
        default = 'NNNXXXXNN',
        help = glue(
            "Several techniques that use UMIs mix the UMI sequence in with a library barcode. In this case we want to remove the random part of the barcodes, but leave the library part so that the reads can be de-multiplexed. We specify this using the `--bcPattern` or short flag `-a` parameter to extract. Ns represent the random part of the barcode and Xs the fixed part. For example, in a standard iCLIP experiment, the barcode is made of 3 random bases, followed by a 4 base library barcode, followed by 2 more random bases. Thus the `--bcPattern` or short flag `-a` would be 'NNNXXXXNN'. [default %default]. If %extractMethod regex, bcPattern must be a sequence of bases as adapter.",  .sep = "\n"
        ),
        dest = "bcPattern"
    ),
    make_option(
        opt_str = c("-s", "--mismatch"),
        type = "integer",
        default = 2,
        help = "How many mismatch are allowed in the adapter (bcPattern)  [default %default]",
        dest = "misMatchAllowed"
    ),
    make_option(
        opt_str = c('-u', '--umiLength'),
        type = 'integer',
        default = 12,
        help = "How many bases for the UMI sequence  [default %default]",
        dest = 'umiLength'
    ),
    make_option(
        opt_str = c("-o", "--output"), 
        type = "character", 
        default = "01-umiProcessed",
        help = "Output folder [default %default]",
        dest = "output"),
    make_option(
        opt_str = c("-p", "--processors"), 
        type = "integer", 
        default = 8,
        help = "Number of processors to use [default %default]",
        dest = "procs"
    ),
    make_option(
        opt_str = c("-q", "--sampleprocs"), 
        type = "integer", 
        default = 2,
        help = "Number of samples to process at each time [default %default]",
        dest = "sampleToprocs"
    ),

    # make_option(
    #     opt_str = c('-w', '--window'), 
    #     type = 'integer', 
    #     default = 10,
    #     help = 'Quality window to use during trimming. Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold. [default %default]',
    #     dest = 'window'),
    # make_option(
    #     opt_str = c('-L', '--leading'), 
    #     type = 'integer', 
    #     default = 3,
    #     help = 'Remove leading low quality or N bases. Cut bases off the start of a read, if below a threshold quality [default %default]',
    #     dest = 'leading'
    # ),
    # make_option(
    #     opt_str = c('-t', '--trailing'), 
    #     type = 'integer', 
    #     default = 3,
    #     help = 'Remove trailing low quality or N bases. Cut bases off the end of a read, if below a threshold quality [default %default]',
    #     dest = 'trailing'
    # ),
    make_option(
        opt_str = c("-z", "--libraryType"),
        type  = 'character', 
        default = "pairEnd",
        help = "The library type to use. Available: 'pairEnd' or 'singleEnd'. [ default %default]",
        dest = "libraryType"
    )#,
    # make_option(
    #     opt_str = c("-m", "--miniumumLength"), 
    #     type = "integer", 
    #     default = 50,
    #     help = "Discard reads less then minimum length [default %default]",
    #     dest = "minL"
    # )
    
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list, description =  paste('Author: OLIVEIRA, H.C.', 'Version: 0.3.5', 'E-mail: hanielcedraz@gmail.com', sep = "\n", collapse = '\n')))






########################################
### CHECKING IF umi_tools EXIST IN THE SYSTEM
########################################

# system("conda activate")
# system("conda install -c bioconda umi_tools")
# pwdInstal <- paste0(getwd(), "/programs/")
# umiToolsInstalled <- Sys.which(paste0(pwdInstal,"UMI-Tools/bin/umi_tools"))
# pipInstalled <- Sys.which("pip3")
# if (opt$installUMItools) {
#     if (nchar(pipInstalled) == 0) {
#         write(paste("pip is not installed. pip is required in order to install UMI-tools. Installing pip..."), stderr())
#         system("python3 -m pip install --upgrade pip")
#     } 
#     if (nchar(umiToolsInstalled) == 0) {
#         write("umi_tools is dependent on python>=3.5, numpy, pandas, scipy, cython, pysam, future, regex and matplotlib", stdout())
#         write(paste("Installing UMItools"), stderr())
#         #system("pip3 install umi_tools --user")
#         system(paste0("pip3 install --target=", pwdInstal, "UMI-Tools/ umi_tools --no-warn-script-location"))
#         cat("/n")
#         write(paste0("umi_tools was installed in ", pwdInstal, "UMI-Tools/bin/"), stderr())
#         umi_tools <- paste(pwdInstal, "UMI-Tools/bin/umi_tools")
#     }
# }


# if (detectCores() < opt$procs) {
#     write(paste("number of cores specified (", opt$procs,") is greater than the number of cores available (",detectCores(),")"), stdout())
#     paste('Using ', detectCores(), 'threads')
# }

# verify if sample_file exist
if (!file.exists(opt$samplesFile)) {
    write(paste("Sample file", opt$samplesFile, "does not exist\n"), stderr())
    stop(paste("Could not find", opt$samplesFile))
}

if (!any(opt$libraryType %in% c('pairEnd', 'singleEnd'))) {
    stop(glue("Library type ({opt$libraryType}) not found. Please provide a valid library type (pairEnd or singleEnd)"))
}

#opt$singleEnd <- TRUE


samples <- baqcomPackage::loadSamplesFile(
    file = opt$samplesFile, 
    reads_folder = opt$rawFolder, 
    column = opt$samplesColumn, libraryType = opt$libraryType
)
cat("####===============================================================\n List of samples\n")
print(samples)
cat("\n")
procs <- baqcomPackage::prepareCore(nThreads = opt$procs)
cat("\n")
#procs <- prepareCore(nThreads = opt$procs)

# create report folders
logFolder <- 'logFolder'
if (!file.exists(file.path(logFolder))) dir.create(file.path(logFolder), recursive = TRUE, showWarnings = FALSE)

# report_folder <- 'reportmiRbaqcomQC-umiTools'
# if (!file.exists(file.path(paste0(reports,'/',report_folder)))) dir.create(file.path(paste0(reports,'/',report_folder)), recursive = TRUE, showWarnings = FALSE)



# 
# #Creating Fastqc plots before Quality Control
# beforeQC <- 'FastQCBefore'
# if (opt$fastqc) {
#     write(paste("Start fastqc - FastQCBefore"), stderr())
#     if (!file.exists(file.path(paste0(reports,'/', beforeQC)))) dir.create(file.path(paste0(reports,'/', beforeQC)), recursive = TRUE, showWarnings = FALSE)
#     fastq.defore <- mclapply(qcquery, function(index){
#         try({
#             system(
#                 paste('fastqc',
#                       paste0(opt$rawFolder, '/',index$sampleName,'*'),
#                       #paste0(opt$rawFolder, '/',index$R2),
#                       ' -o ',
#                       paste0(reports,'/', beforeQC),
#                       '-t', opt$sampleToprocs
#                 )
#             )
#         })
#     }, mc.cores = opt$sampleToprocs)
#     # for (i in samples[,1]) {
#     #         system2('fastqc',
#     #                paste0(opt$rawFolder, '/', i,'*', ' -o ', paste0(reports,'/', beforeQC), ' -t ', opt$sampleToprocs))
#     #     }
# }






## create output folder
output_Folder <- opt$output
if (opt$umiCommand == "extract") {
    extractedFolder <- paste0(output_Folder, "/00-extracted")
    if (!file.exists(file.path(extractedFolder))) {
        dir.create(file.path(extractedFolder), recursive = TRUE, showWarnings = FALSE)
    }
    write(glue("Creating directory: {extractedFolder}"), stdout())
} else if (opt$umiCommand == "dedup") {
    dedupFolder <- paste0(output_Folder, "/01-deduplicated")
    if (!file.exists(file.path(dedupFolder))) {
        dir.create(file.path(dedupFolder), recursive = TRUE, showWarnings = FALSE)
    }
    write(glue("Creating directory: {dedupFolder}"), stdout())
} else if (opt$umiCommand == "both") {
    extractedFolder <- paste0(output_Folder, "/00-extracted")
    if (!file.exists(file.path(extractedFolder))) {
        dir.create(file.path(extractedFolder), recursive = TRUE, showWarnings = FALSE)
    }
    bowtieFolder <- paste0(output_Folder, "/01-bowtieMapped")
    if (!file.exists(file.path(bowtieFolder))) {
        dir.create(file.path(bowtieFolder), recursive = TRUE, showWarnings = FALSE)
    }
    dedupFolder <- paste0(output_Folder, "/02-deduplicated")
    if (!file.exists(file.path(dedupFolder))) {
        dir.create(file.path(dedupFolder), recursive = TRUE, showWarnings = FALSE)
    }
    write(glue("Creating directory: {extractedFolder}."), stdout())
    write(glue("Creating directory: {bowtieFolder}."), stdout())
    write(glue("Creating directory: {dedupFolder}."), stdout())
}




# UMI-Tools analysis function
#pigz <- system('which pigz 2> /dev/null', ignore.stdout = TRUE, ignore.stderr = TRUE)

if (opt$umiCommand == "extract") {
    umiQuery <- baqcomPackage::createSampleList(
        samples = samples, 
        reads_folder = opt$rawFolder, 
        column = opt$samplesColumn, fileType = "fastq.gz", 
        libraryType = opt$libraryType, 
        program = "trimmomatic"
    )
    opt$extractMethod <- "regex"
    opt$discard <- FALSE
    opt$bcPattern <- "AACTGTAGGCACCATCAAT"
    if (opt$extractMethod == "regex") {
        if (opt$discard) {
            opt$bcPattern <- paste0("'.+(?P<discard_1>",opt$bcPattern, "){s<=", opt$misMatchAllowed, "}(?P<umi_1>.{",opt$umiLength,"})'")
        } else {
            opt$bcPattern <- paste0("'.+(",opt$bcPattern, "){s<=", opt$misMatchAllowed, "}(?P<umi_1>.{",opt$umiLength,"})'")
        }
        
    }
    
    cat("####===============================================================\n List of samples to process\n")
    print(umiQuery)
    cat("\n")
    if (opt$libraryType == "pairEnd") {
        cat("\nStarting umi extract for Peired-End library\n\n")
        umiPair <- mclapply(umiQuery, function(index) {
            try({
                system(
                    paste(
                        "umi_tools",
                        "extract",
                        paste0("-I", paste0(opt$rawFolder, "/", index$R1)),
                        paste0("--extract-method=", opt$extractMethod),
                        paste0("--bc-pattern", opt$bcPattern, sep = "="),
                        paste0("--read2-in", paste0(opt$rawFolder, "/", index$R2), sep = "="),
                        paste0("--stdout", paste0(extractedFolder, "/", index$sampleName, "_extracted_PE1.fastq.gz"), sep = "="),
                        paste("--read2-out", paste0(extractedFolder, "/", index$sampleName, "_extracted_PE2.fastq.gz"), sep = "=")
                    )
                )
            })
        }, mc.cores = opt$sampleToprocs)
        
        if (!all(sapply(umiPair , "==", 0L))) {
            write(paste("Something went wrong with umi_tools some jobs failed"),stderr())
            stop()
        }
        
        write(glue("Results from umi extract were stored in {extractedFolder}"), stdout())
        cat("\n=============================\n")
        write(glue("Log files from umi extract were stored in {extractedFolder}"), stdout())
        cat("\n=============================\n")
        
    } else if (opt$libraryType == "singleEnd") {
        cat("\nStarting umi extract for Single-End library\n\n")
        
        umisingle <- mclapply(umiQuery, function(index) {
            write(paste(
                "umi_tools",
                "extract",
                paste0("--stdin=", index$SE),
                paste0("--extract-method=", opt$extractMethod),
                paste0("--bc-pattern=", opt$bcPattern),
                paste0("--log=", paste0(logFolder, "/", index$sampleName, "_extracted.log")),
                paste("--stdout=", paste0(extractedFolder, "/", index$sampleName, "_extracted_SE.fastq.gz"))
            ), stdout())
            cat("\n")
            
            try({
                system(
                    paste(
                        "umi_tools",
                        "extract",
                        paste0("--stdin=", index$SE),
                        paste0("--extract-method=", opt$extractMethod),
                        paste0("--bc-pattern=", opt$bcPattern),
                        paste0("--log=", paste0(logFolder, "/", index$sampleName, "_extracted.log")),
                        paste("--stdout=", paste0(extractedFolder, "/", index$sampleName, "_extracted_SE.fastq.gz"))
                    )
                )
            })
        }, mc.cores = opt$sampleToprocs)
        
        if (!all(sapply(umisingle , "==", 0L))) {
            write(paste("Something went wrong with umi_tools some jobs failed"),stderr())
            stop()
        }
        
        write(glue("Results from umi extract were stored in {extractedFolder}"), stdout())
        cat("\n=============================\n")
        write(glue("Log files from umi extract were stored in {extractedFolder}"), stdout())
        cat("\n=============================\n")
    }
} else if (opt$umiCommand == "dedup") {
    umiQuery <- baqcomPackage::createSampleList(
        samples = samples, 
        reads_folder = opt$rawFolder, 
        column = opt$samplesColumn, fileType = "bam", 
        libraryType = opt$libraryType, 
        program = "trimmomatic"
    )

    cat("####===============================================================\n List of samples to process\n")
    print(umiQuery)
    cat("\n")
    umidedup <- mclapply(umiQuery, function(index) {
        try({
            system(
                paste(
                    "umi_tools",
                    "dedup",
                    paste("--stdin", paste0(opt$rawFolder, "/", index$SE)),
                    paste("--bc-pattern", opt$bcPattern, sep = "="),
                    paste("--log", paste0(logFolder, "/", index$sampleName, "_extracted.log"), sep = "="),
                    paste("--stdout", paste0(extractedFolder, "/", index$sampleName, "_extracted_SE.fastq.gz"), sep = "=")
                )
            )
        })
    }, mc.cores = opt$sampleToprocs)
    
    if (!all(sapply(umidedup , "==", 0L))) {
        write(paste("Something went wrong with umi_tools some jobs failed"),stderr())
        stop()
    }
    
} 




#umi_tools dedup -I example.bam --output-stats=deduplicated -S deduplicated.bam

# 
# 
# if (!opt$singleEnd) {
#     write(paste("START Trimmomatic PE"), stderr())
#     trimmomatic.pair <- mclapply(umiQuery, function(index){
#         
#         try({
#             system(
#                 paste('java', '-jar', trimmomatic, 'PE',
#                       paste0('-threads ',
#                              ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs))),
#                       paste0(opt$rawFolder, '/', index$R1),
#                       paste0(opt$rawFolder, '/', index$R2),
#                       if (pigz == 0) {
#                           paste(paste0(opt$output, '/', index$sampleName, '_', 'trim_PE1.fastq'),
#                                 paste0(opt$output, '/', index$sampleName, '_', 'trim_SE1.fastq'),
#                                 paste0(opt$output, '/', index$sampleName, '_', 'trim_PE2.fastq'),
#                                 paste0(opt$output, '/', index$sampleName, '_', 'trim_SE2.fastq'))}
#                       else{
#                           paste(paste0(opt$output, '/', index$sampleName, '_', 'trim_PE1.fastq.gz'),
#                                 paste0(opt$output, '/', index$sampleName, '_', 'trim_SE1.fastq.gz'),
#                                 paste0(opt$output, '/', index$sampleName, '_', 'trim_PE2.fastq.gz'),
#                                 paste0(opt$output, '/', index$sampleName, '_', 'trim_SE2.fastq.gz'))},
#                       paste0('-summary ', reports,'/', report_folder, '/', index$sampleName, '_', 'statsSummaryFile.txt'),
#                       paste0('ILLUMINACLIP:', trimmomatic_dir, 'adapters', '/',
#                              opt$adapters, ':2:30:10'),
#                       paste0('LEADING:', opt$leading),
#                       paste0('TRAILING:', opt$trailing),
#                       paste0('SLIDINGWINDOW:', opt$window, ':', opt$qual),
#                       paste0('MINLEN:', opt$minL),
#                       paste0('2> ', reports,'/',report_folder, '/', index$sampleName, '_', 'trimmomatic_out.log'), collapse = '\t'), ignore.stdout = TRUE
#             )
#         })
#     }, mc.cores = opt$sampleToprocs)
#     
#     if (!all(sapply(trimmomatic.pair , "==", 0L))) {
#         write(paste("Something went wrong with Trimmomatic some jobs failed"),stderr())
#         stop()
#     }
# } else if (opt$singleEnd) {
#     write(paste("START Trimmomatic SE"), stderr())
#     trimmomatic.single <- mclapply(qcquery, function(index){
#         try({
#             system(
#                 paste('java', '-jar', trimmomatic, 'SE',
#                       paste0('-threads ',
#                              ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs))),
#                       paste0(opt$rawFolder, '/', index$SE),
#                       if (pigz == 0) {
#                           paste0(opt$output, '/', index$sampleName, '_', 'trim_SE.fastq')
#                       }
#                       else{
#                           paste0(opt$output, '/', index$sampleName, '_', 'trim_SE.fastq.gz')
#                       },
#                       paste0('-summary ', reports,'/', report_folder, '/', index$sampleName, '_', 'statsSummaryFile.txt'),
#                       paste0('ILLUMINACLIP:', trimmomatic_dir, 'adapters', '/',
#                              opt$adapters, ':2:30:10'),
#                       paste0('LEADING:', opt$leading),
#                       paste0('TRAILING:', opt$trailing),
#                       paste0('SLIDINGWINDOW:', opt$window, ':', opt$qual),
#                       paste0('MINLEN:', opt$minL),
#                       paste0('2> ', reports,'/',report_folder, '/', index$sampleName, '_', 'trimmomatic_out.log'), collapse = '\t'), ignore.stdout = TRUE
#             )
#         })
#     }, mc.cores = opt$sampleToprocs)
#     
#     # if (!all(sapply(trimmomatic.single , "==", 0L))) {
#     #     write(paste("Something went wrong with Trimmomatic some jobs failed"),stderr())
#     #     stop()
#     # }
# }
# 
# if (pigz == 0) {
#     write(paste("Compressing files with pigz using", procs, 'processors'),stderr())
#     system(paste0('pigz ', '-f ', '-p ', procs,' ', opt$output,'/*.fastq'))
# }
# 
# 
# cat('\n')

# output_folder <- opt$output
# if (!opt$singleEnd) {
#     afterqcList <- function(samples, output_folder, column){
#         mapping_list <- list()
#         for (i in 1:nrow(samples)) {
#             reads <- dir(path = file.path(output_folder), pattern = "fastq.gz$", full.names = TRUE)
#             # for (i in seq.int(to=nrow(samples))){
#             #     reads <- dir(path=file.path(reads_folder,samples[i,column]),pattern="gz$",full.names=TRUE)
#             map <- lapply(c("_PE1", "_PE2", "_SE1", "_SE2"), grep, x = reads, value = TRUE)
#             names(map) <- c("PE1", "PE2", "SE1", "SE2")
#             map$sampleName <-  samples[i,column]
#             map$PE1 <- map$PE1[i]
#             map$PE2 <- map$PE2[i]
#             map$SE1 <- map$SE1[i]
#             map$SE2 <- map$SE2[i]
#             mapping_list[[paste(map$sampleName)]] <- map
#             mapping_list[[paste(map$sampleName, sep = "_")]]
#         }
#         write(paste("Setting up", length(mapping_list), "jobs"),stdout())
#         return(mapping_list)
#     }
#     
# } else if (opt$singleEnd) {
#     afterqcList <- function(samples, reads_folder, column){
#         mapping_list <- list()
#         for (i in 1:nrow(samples)) {
#             reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
#             #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
#             map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
#             names(map) <- c("SE")
#             map$sampleName <-  samples[i,column]
#             map$SE <- map$SE[i]
#             #map$R2 <- samples[i,3]
#             mapping_list[[paste(map$sampleName)]] <- map
#             #mapping_list[[paste(map$sampleName)]]
#         }
#         write(paste("Setting up",length(mapping_list),"jobs"), stdout())
#         return(mapping_list)
#     }
# }
# 
# query.after <- afterqcList(samples, opt$output, opt$samplesColumn)
# 
# #Creating Fastqc plots after Quality Control
# afterQC <- 'FastQCAfter'
# if (opt$fastqc) {
#     write(paste("Start fastqc - FastQCAfter"), stderr())
#     if (!file.exists(file.path(paste0(reports,'/', afterQC)))) dir.create(file.path(paste0(reports,'/', afterQC)), recursive = TRUE, showWarnings = FALSE)
#     fastq.after <- mclapply(qcquery, function(index){
#         try({
#             system(
#                 paste('fastqc',
#                       paste0(opt$output, '/',index$sampleName,'*'),
#                       #paste0(opt$rawFolder, '/',index$R2),
#                       ' -o ',
#                       paste0(reports,'/', afterQC),
#                       '-t', opt$sampleToprocs
#                 )
#             )
#         })
#     }, mc.cores = opt$sampleToprocs)
#     # for (i in samples[,1]) {
#     #     system2('fastqc',
#     #             paste0(opt$output, '/', i, '_trim_PE1*', ' -o ', paste0(reports,'/', afterQC), ' -t ', opt$sampleToprocs))
#     #     system2('fastqc',
#     #             paste0(opt$output, '/', i, '_trim_PE2*', ' -o ', paste0(reports,'/',afterQC), ' -t ', opt$sampleToprocs))
#     #     }
# }
# 
# cat('\n')
# 
# if (opt$multiqc) {
#     if (file.exists(paste0(reports,'/', beforeQC)) & file.exists(paste0(reports,'/', afterQC))) {
#         system2('multiqc', paste(paste0(reports,'/', report_folder), paste0(reports,'/',beforeQC), paste0(reports,'/',afterQC), '-o', reports, '-f'))
#     }else{
#         system2('multiqc', paste(paste0(reports,'/', report_folder), '-o', reports, '-f'))
#     }
# }

# cat('\n')
# 
# # Creating samples report
# if (!opt$singleEnd) {
#     TidyTable <- function(x) {
#         final <- data.frame('Input_Read_Pairs' = x[1,2],
#                             'Pairs_Reads' = x[2,2],
#                             'Pairs_Read_Percent' = x[3,2],
#                             'Forward_Only_Surviving_Reads' = x[4,2],
#                             'Forward_Only_Surviving_Read_Percent' = x[5,2],
#                             'Reverse_Only_Surviving_Reads' = x[6,2],
#                             'Reverse_Only_Surviving_Read_Percent' = x[7,2],
#                             'Dropped_Reads' = x[8,2],
#                             'Dropped_Read_Percent' = x[9,2])
#         return(final)
#     }
#     
#     report_sample <- list()
#     for (i in samples[,1]) { # change this to your "samples"
#         report_sample[[i]] <- read.table(paste0(reports,'/', report_folder, '/', i,"_statsSummaryFile.txt"),
#                                          header = F, as.is = T, fill = TRUE, sep = ':', text = TRUE)
#     }
#     
#     df <- lapply(report_sample, FUN = function(x) TidyTable(x))
#     final_df <- do.call("rbind", df)
# }else if (opt$singleEnd) {
#     TidyTable <- function(x) {
#         final <- data.frame('Input_Read_Pairs' = x[1,2],
#                             # 'Pairs_Reads' = x[2,2],
#                             # 'Pairs_Read_Percent' = x[3,2],
#                             'Surviving_Reads' = x[2,2],
#                             'Surviving_Read_Percent' = x[3,2],
#                             # 'Reverse_Only_Surviving_Reads' = x[6,2],
#                             # 'Reverse_Only_Surviving_Read_Percent' = x[7,2],
#                             'Dropped_Reads' = x[4,2],
#                             'Dropped_Read_Percent' = x[5,2]
#         )
#         return(final)
#     }
#     
#     report_sample <- list()
#     for (i in samples[,1]) { # change this to your "samples"
#         report_sample[[i]] <- read.table(paste0(reports,'/', report_folder, '/', i,"_statsSummaryFile.txt"),
#                                          header = F, as.is = T, fill = TRUE, sep = ':', text = TRUE)
#     }
#     
#     df <- lapply(report_sample, FUN = function(x) TidyTable(x))
#     final_df <- do.call("rbind", df)
# }
# 
# write.table(final_df, file = paste0(reports, '/', 'QualityControlReportSummary.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)
# 
# system2('cat', paste0(reports,'/','QualityControlReportSummary.txt'))
# 
cat('\n')
write(paste('How to cite:', sep = '\n', collapse = '\n', "Please, visit https://github.com/hanielcedraz/BAQCOM/blob/master/how_to_cite.txt", "or see the file 'how_to_cite.txt'"), stderr())

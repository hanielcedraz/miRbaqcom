#!/usr/bin/env Rscript


########################################
### LOADING PACKAGES
########################################
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("glue"))

########################################
### SETING PARAMETERS
########################################
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")


option_list <- list(
    make_option(
        opt_str = c("-f", "--file"), 
        type = "character", 
        default = "samples.txt",
        help = "The filename of the sample file [default %default]",
        dest = "samplesFile"
    ),
    make_option(
        opt_str = c("-d", "--directory"), 
        type = "character", 
        default = "00-Fastq",
        help = "Directory where the raw sequence data is stored [default %default]",
        dest = "Raw_Folder"
    ),
    make_option(
        opt_str = c("-c", "--column"), 
        type = "character", 
        default = "SAMPLE_ID",
        help = "Column name from the sample sheet to use as read folder names [default %default]",
        dest = "samplesColumn"
    ),
    make_option(
        opt_str = c("-a", "--adapters"), 
        type  = 'character', 
        default = 'AACCGGTT',
        help = glue(
            "Sequences of adapters ligated to the 3' and 5' end (paired data: of the first read).", "For the 3' end, the adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found if it is a suffix of the read. [default %default],",  .sep = "\n"
        ),
        dest = "adapters"
    ),
    make_option(
        opt_str = c("-r", "--multiqc"), 
        action = 'store_true', 
        type = "logical",
        default = FALSE,
        help = "Use this option if you want to run multiqc software  [default %default]",
        dest = "multiqc"
    ),
    make_option(
        opt_str = c("-o", "--output"), 
        type = "character", 
        default = "01-CleanedReads",
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
    make_option(
        opt_str = c('-s', '--sliding'), 
        type = 'integer', 
        default = 20,
        help = 'Quality sliding to use during trimming. Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold. [default %default]',
        dest = 'qual'
    ),
    make_option(
        opt_str = c('-w', '--window'), 
        type = 'integer', 
        default = 10,
        help = 'Quality window to use during trimming. Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold. [default %default]',
        dest = 'window'),
    make_option(
        opt_str = c('-L', '--leading'), 
        type = 'integer', 
        default = 3,
        help = 'Remove leading low quality or N bases. Cut bases off the start of a read, if below a threshold quality [default %default]',
        dest = 'leading'
    ),
    make_option(
        opt_str = c('-t', '--trailing'), 
        type = 'integer', 
        default = 3,
        help = 'Remove trailing low quality or N bases. Cut bases off the end of a read, if below a threshold quality [default %default]',
        dest = 'trailing'
    ),
    make_option(
        opt_str = c("-z", "--single"), 
        action = "store_true", 
        default = FALSE,
        help = "Use this option if you have single-end files [doesn't need an argument]. [%default]",
        dest = "singleEnd"
    ),
    make_option(
        opt_str = c("-m", "--miniumumLength"), 
        type = "integer", 
        default = 50,
        help = "Discard reads less then minimum length [default %default]",
        dest = "minL"
    )
    
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list, description =  paste('Author: OLIVEIRA, H.C.', 'Version: 0.3.5', 'E-mail: hanielcedraz@gmail.com', sep = "\n", collapse = '\n')))



########################################
### PERFORMING QC ANALYSIS
########################################


multiqc <- system('which multiqc 2> /dev/null', ignore.stdout = TRUE, ignore.stderr = TRUE)
if (opt$multiqc) {
    if (multiqc != 0) {
        write(paste("Multiqc is not installed. If you would like to use multiqc analysis, please install it or remove -r parameter"), stderr())
        stop()
    }
}


if (detectCores() < opt$procs) {
    write(paste("number of cores specified (", opt$procs,") is greater than the number of cores available (",detectCores(),")"), stdout())
    paste('Using ', detectCores(), 'threads')
}

# verify if sample_file exist
if ( !file.exists(opt$samplesFile) ) {
    write(paste("Sample file",opt$samplesFile,"does not exist\n"), stderr())
    stop()
}


## loadSampleFile
loadSamplesFile <- function(file, reads_folder, column){
    ## debug
    file = opt$samplesFile; reads_folder = opt$Raw_Folder; column = opt$samplesColumn
    ##
    if ( !file.exists(file) ) {
        write(paste("Sample file",file,"does not exist\n"), stderr())
        stop()
    }
    ### column SAMPLE_ID should be the sample name
    ### rows can be commented out with #
    targets <- read.table(file, sep = "", header = TRUE, as.is = TRUE)
    if (!opt$singleEnd) {
        if (!all(c("SAMPLE_ID", "Read_1", "Read_2") %in% colnames(targets))) {
            write(paste("Expecting the three columns SAMPLE_ID, Read_1 and Read_2 in samples file (tab-delimited)\n"), stderr())
            stop()
        }
    }
    for (i in seq.int(nrow(targets$SAMPLE_ID))) {
        if (targets[i, column]) {
            ext <- unique(file_ext(dir(file.path(reads_folder, targets[i,column]), pattern = "gz")))
            if (length(ext) == 0) {
                write(paste("Cannot locate fastq or sff file in folder",targets[i,column],"\n"), stderr())
                stop()
            }
            # targets$type[i] <- paste(ext,sep="/")
        }
        else {
            ext <- file_ext(grep("gz", dir(file.path(reads_folder, targets[i, column])), value = TRUE))
            if (length(ext) == 0) {
                write(paste(targets[i,column],"is not a gz file\n"), stderr())
                stop()
            }
            
        }
    }
    write(paste("samples.txt contains", nrow(targets), "samples to process"),stdout())
    return(targets)
}


prepareCore <- function(opt_procs){
    # if opt_procs set to 0 then expand to samples by targets
    if (detectCores() < opt$procs) opt_procs <- detectCores()
    write(paste("Using",opt_procs, "processors", sep = " "),stdout())
    return(opt_procs)
}

if (!opt$singleEnd) {
    qcList <- function(samples, reads_folder, column){
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
    }
} else if (opt$singleEnd) {
    qcList <- function(samples, reads_folder, column){
        mapping_list <- list()
        for (i in 1:nrow(samples)) {
            reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
            #reads <- dir(path=file.path(reads_folder, samples[i,column]), pattern = "fastq.gz$", full.names = TRUE)
            map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
            names(map) <- c("SE")
            map$sampleName <-  samples[i,column]
            #map$SE <- map$SE[i]
            map$SE <- samples[i,2]
            mapping_list[[paste(map$sampleName)]] <- map
            #mapping_list[[paste(map$sampleName)]]
        }
        write(paste("Setting up",length(mapping_list),"jobs"), stdout())
        return(mapping_list)
    }
}

samples <- loadSamplesFile(opt$samplesFile, opt$Raw_Folder, opt$samplesColumn)
procs <- prepareCore(opt$procs)
qcquery <- qcList(samples, opt$Raw_Folder, opt$samplesColumn)


# create report folders
reports <- '02-Reports'
if (!file.exists(file.path(reports))) dir.create(file.path(reports), recursive = TRUE, showWarnings = FALSE)

report_folder <- 'reportBaqcomQC'
if (!file.exists(file.path(paste0(reports,'/',report_folder)))) dir.create(file.path(paste0(reports,'/',report_folder)), recursive = TRUE, showWarnings = FALSE)



#Creating Fastqc plots before Quality Control
beforeQC <- 'FastQCBefore'
if (opt$fastqc) {
    write(paste("Start fastqc - FastQCBefore"), stderr())
    if (!file.exists(file.path(paste0(reports,'/', beforeQC)))) dir.create(file.path(paste0(reports,'/', beforeQC)), recursive = TRUE, showWarnings = FALSE)
    fastq.defore <- mclapply(qcquery, function(index){
        try({
            system(
                paste('fastqc',
                      paste0(opt$Raw_Folder, '/',index$sampleName,'*'),
                      #paste0(opt$Raw_Folder, '/',index$R2),
                      ' -o ',
                      paste0(reports,'/', beforeQC),
                      '-t', opt$sampleToprocs
                )
            )
        })
    }, mc.cores = opt$sampleToprocs)
    # for (i in samples[,1]) {
    #         system2('fastqc',
    #                paste0(opt$Raw_Folder, '/', i,'*', ' -o ', paste0(reports,'/', beforeQC), ' -t ', opt$sampleToprocs))
    #     }
}






## create output folder
output_Folder <- opt$output
if (!file.exists(file.path(output_Folder))) dir.create(file.path(output_Folder), recursive = TRUE, showWarnings = FALSE)



# Trimmomatic analysis function
pigz <- system('which pigz 2> /dev/null', ignore.stdout = TRUE, ignore.stderr = TRUE)
if (!opt$singleEnd) {
    write(paste("START Trimmomatic PE"), stderr())
    trimmomatic.pair <- mclapply(qcquery, function(index){
        try({
            system(
                paste('java', '-jar', trimmomatic, 'PE',
                      paste0('-threads ',
                             ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs))),
                      paste0(opt$Raw_Folder, '/', index$R1),
                      paste0(opt$Raw_Folder, '/', index$R2),
                      if (pigz == 0) {
                          paste(paste0(opt$output, '/', index$sampleName, '_', 'trim_PE1.fastq'),
                                paste0(opt$output, '/', index$sampleName, '_', 'trim_SE1.fastq'),
                                paste0(opt$output, '/', index$sampleName, '_', 'trim_PE2.fastq'),
                                paste0(opt$output, '/', index$sampleName, '_', 'trim_SE2.fastq'))}
                      else{
                          paste(paste0(opt$output, '/', index$sampleName, '_', 'trim_PE1.fastq.gz'),
                                paste0(opt$output, '/', index$sampleName, '_', 'trim_SE1.fastq.gz'),
                                paste0(opt$output, '/', index$sampleName, '_', 'trim_PE2.fastq.gz'),
                                paste0(opt$output, '/', index$sampleName, '_', 'trim_SE2.fastq.gz'))},
                      paste0('-summary ', reports,'/', report_folder, '/', index$sampleName, '_', 'statsSummaryFile.txt'),
                      paste0('ILLUMINACLIP:', trimmomatic_dir, 'adapters', '/',
                             opt$adapters, ':2:30:10'),
                      paste0('LEADING:', opt$leading),
                      paste0('TRAILING:', opt$trailing),
                      paste0('SLIDINGWINDOW:', opt$window, ':', opt$qual),
                      paste0('MINLEN:', opt$minL),
                      paste0('2> ', reports,'/',report_folder, '/', index$sampleName, '_', 'trimmomatic_out.log'), collapse = '\t'), ignore.stdout = TRUE
            )
        })
    }, mc.cores = opt$sampleToprocs)
    
    if (!all(sapply(trimmomatic.pair , "==", 0L))) {
        write(paste("Something went wrong with Trimmomatic some jobs failed"),stderr())
        stop()
    }
} else if (opt$singleEnd) {
    write(paste("START Trimmomatic SE"), stderr())
    trimmomatic.single <- mclapply(qcquery, function(index){
        try({
            system(
                paste('java', '-jar', trimmomatic, 'SE',
                      paste0('-threads ',
                             ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs))),
                      paste0(opt$Raw_Folder, '/', index$SE),
                      if (pigz == 0) {
                          paste0(opt$output, '/', index$sampleName, '_', 'trim_SE.fastq')
                      }
                      else{
                          paste0(opt$output, '/', index$sampleName, '_', 'trim_SE.fastq.gz')
                      },
                      paste0('-summary ', reports,'/', report_folder, '/', index$sampleName, '_', 'statsSummaryFile.txt'),
                      paste0('ILLUMINACLIP:', trimmomatic_dir, 'adapters', '/',
                             opt$adapters, ':2:30:10'),
                      paste0('LEADING:', opt$leading),
                      paste0('TRAILING:', opt$trailing),
                      paste0('SLIDINGWINDOW:', opt$window, ':', opt$qual),
                      paste0('MINLEN:', opt$minL),
                      paste0('2> ', reports,'/',report_folder, '/', index$sampleName, '_', 'trimmomatic_out.log'), collapse = '\t'), ignore.stdout = TRUE
            )
        })
    }, mc.cores = opt$sampleToprocs)
    
    # if (!all(sapply(trimmomatic.single , "==", 0L))) {
    #     write(paste("Something went wrong with Trimmomatic some jobs failed"),stderr())
    #     stop()
    # }
}

if (pigz == 0) {
    write(paste("Compressing files with pigz using", procs, 'processors'),stderr())
    system(paste0('pigz ', '-f ', '-p ', procs,' ', opt$output,'/*.fastq'))
}


cat('\n')

output_folder <- opt$output
if (!opt$singleEnd) {
    afterqcList <- function(samples, output_folder, column){
        mapping_list <- list()
        for (i in 1:nrow(samples)) {
            reads <- dir(path = file.path(output_folder), pattern = "fastq.gz$", full.names = TRUE)
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
    }
    
} else if (opt$singleEnd) {
    afterqcList <- function(samples, reads_folder, column){
        mapping_list <- list()
        for (i in 1:nrow(samples)) {
            reads <- dir(path = file.path(reads_folder), pattern = "fastq.gz$", full.names = TRUE)
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
    }
}

query.after <- afterqcList(samples, opt$output, opt$samplesColumn)

#Creating Fastqc plots after Quality Control
afterQC <- 'FastQCAfter'
if (opt$fastqc) {
    write(paste("Start fastqc - FastQCAfter"), stderr())
    if (!file.exists(file.path(paste0(reports,'/', afterQC)))) dir.create(file.path(paste0(reports,'/', afterQC)), recursive = TRUE, showWarnings = FALSE)
    fastq.after <- mclapply(qcquery, function(index){
        try({
            system(
                paste('fastqc',
                      paste0(opt$output, '/',index$sampleName,'*'),
                      #paste0(opt$Raw_Folder, '/',index$R2),
                      ' -o ',
                      paste0(reports,'/', afterQC),
                      '-t', opt$sampleToprocs
                )
            )
        })
    }, mc.cores = opt$sampleToprocs)
    # for (i in samples[,1]) {
    #     system2('fastqc',
    #             paste0(opt$output, '/', i, '_trim_PE1*', ' -o ', paste0(reports,'/', afterQC), ' -t ', opt$sampleToprocs))
    #     system2('fastqc',
    #             paste0(opt$output, '/', i, '_trim_PE2*', ' -o ', paste0(reports,'/',afterQC), ' -t ', opt$sampleToprocs))
    #     }
}

cat('\n')

if (opt$multiqc) {
    if (file.exists(paste0(reports,'/', beforeQC)) & file.exists(paste0(reports,'/', afterQC))) {
        system2('multiqc', paste(paste0(reports,'/', report_folder), paste0(reports,'/',beforeQC), paste0(reports,'/',afterQC), '-o', reports, '-f'))
    }else{
        system2('multiqc', paste(paste0(reports,'/', report_folder), '-o', reports, '-f'))
    }
}

cat('\n')

# Creating samples report
if (!opt$singleEnd) {
    TidyTable <- function(x) {
        final <- data.frame('Input_Read_Pairs' = x[1,2],
                            'Pairs_Reads' = x[2,2],
                            'Pairs_Read_Percent' = x[3,2],
                            'Forward_Only_Surviving_Reads' = x[4,2],
                            'Forward_Only_Surviving_Read_Percent' = x[5,2],
                            'Reverse_Only_Surviving_Reads' = x[6,2],
                            'Reverse_Only_Surviving_Read_Percent' = x[7,2],
                            'Dropped_Reads' = x[8,2],
                            'Dropped_Read_Percent' = x[9,2])
        return(final)
    }
    
    report_sample <- list()
    for (i in samples[,1]) { # change this to your "samples"
        report_sample[[i]] <- read.table(paste0(reports,'/', report_folder, '/', i,"_statsSummaryFile.txt"),
                                         header = F, as.is = T, fill = TRUE, sep = ':', text = TRUE)
    }
    
    df <- lapply(report_sample, FUN = function(x) TidyTable(x))
    final_df <- do.call("rbind", df)
}else if (opt$singleEnd) {
    TidyTable <- function(x) {
        final <- data.frame('Input_Read_Pairs' = x[1,2],
                            # 'Pairs_Reads' = x[2,2],
                            # 'Pairs_Read_Percent' = x[3,2],
                            'Surviving_Reads' = x[2,2],
                            'Surviving_Read_Percent' = x[3,2],
                            # 'Reverse_Only_Surviving_Reads' = x[6,2],
                            # 'Reverse_Only_Surviving_Read_Percent' = x[7,2],
                            'Dropped_Reads' = x[4,2],
                            'Dropped_Read_Percent' = x[5,2]
        )
        return(final)
    }
    
    report_sample <- list()
    for (i in samples[,1]) { # change this to your "samples"
        report_sample[[i]] <- read.table(paste0(reports,'/', report_folder, '/', i,"_statsSummaryFile.txt"),
                                         header = F, as.is = T, fill = TRUE, sep = ':', text = TRUE)
    }
    
    df <- lapply(report_sample, FUN = function(x) TidyTable(x))
    final_df <- do.call("rbind", df)
}

write.table(final_df, file = paste0(reports, '/', 'QualityControlReportSummary.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)

system2('cat', paste0(reports,'/','QualityControlReportSummary.txt'))

cat('\n')
write(paste('How to cite:', sep = '\n', collapse = '\n', "Please, visit https://github.com/hanielcedraz/BAQCOM/blob/master/how_to_cite.txt", "or see the file 'how_to_cite.txt'"), stderr())

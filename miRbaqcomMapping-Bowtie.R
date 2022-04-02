#!/usr/bin/env Rscript


# Create a new conda environment with all the r-essentials conda packages built from CRAN:
#     
#     conda create -n r_env r-essentials r-base
# 
# Activate the environment:
#     
#     conda activate r_env
# 
# List the packages in the environment:
#     
#     conda list



# Sys.which(paste0("conda"))
# # system(
# #     paste("conda create -n r_env r-essentials r-base")
# # )
# system(
#     paste("conda init")
# )
# system(
#     paste("conda activate r_env")
# )
# system(
#     paste("conda install -c bioconda bowtie")
# )

if (!require(baqcomPackage, quietly = TRUE)) {
    devtools::install_github(repo = "git@github.com:hanielcedraz/baqcomPackage.git", upgrade = "always", quiet = TRUE, force = TRUE)
}
#suppressPackageStartupMessages(devtools::install_github(repo = "git@github.com:hanielcedraz/baqcomPackage.git", upgrade = "never", quiet = TRUE, force = TRUE))
########################################
### LOADING PACKAGES
########################################
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("baqcomPackage"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))

#source("~/Documents/baqcomPackage/R/createSampleListFunction.R")
########################################
### SETING PARAMETERS
########################################
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")

option_list <- list(
    # make_option(
    #     opt_str = c("-u", "--useConda"),
    #     action = 'store_true', 
    #     type = "logical",
    #     default = FALSE,
    #     help = "Wether to use conda enviroment",
    #     dest = "useConda"
    # ),
    make_option(
        opt_str = c("-f", "--file"), 
        type = "character", 
        default = "samples.txt",
        help = "The filename of the sample file [default %default]",
        dest = "samplesFile"
    ),
    make_option(
        opt_str = c("-d", "--inputDirectory"), 
        type = "character", 
        default = "01-CleanedReads",
        help = "Directory where the raw sequence data is stored [default %default]",
        dest = "cleanedFolder"
    ),
    make_option(
        opt_str = c("-c", "--column"), 
        type = "character", 
        default = "SAMPLE_ID",
        help = "Column name from the sample sheet to use as read folder names [default %default]",
        dest = "samplesColumn"
    ),
    make_option(
        opt_str = c("-m", "--multiqc"), 
        action = 'store_true', 
        type = "logical",
        default = FALSE,
        help = "Use this option if you want to run multiqc software  [default %default]",
        dest = "multiqc"
    ),
    make_option(
        opt_str = c("-o", "--output"), 
        action = "store_true",
        type = "character", 
        #default = '02-MappedReads',
        help = "Directory to store the mapped results [default %default]",
        dest = "mappingFolder"
    ),
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
        opt_str = c("-z", "--libraryType"),
        type  = 'character', 
        default = "singleEnd",
        help = "The library type to use. Available: 'pairEnd' or 'singleEnd'. [ default %default]",
        dest = "libraryType"
    ),
    make_option(
        opt_str = c("-M", "--program"),
        type  = 'character', 
        default = "bowtie",
        help = "Which mapping program to use. Options: 'bowtie', 'bwa'. [ default %default]",
        dest = "mappingProgram"
    ),
    make_option(
        opt_str = c("-t", "--mappingTargets"), 
        type = "character", 
        default = "referenceGenome.fa",
        help = "Path to the fasta file [target fasta] to run mapping against (default %default); or path to the directory where the genome indices are stored (path/to/the/genoma_file/index.",
        dest = "mappingTarget"
    ),
    make_option(
        opt_str = c("-i", "--index"), 
        action = "store_true", 
        default = FALSE,
        help = "This option directs to re-run genome indices generation. [%default]",
        dest = "indexBuild"
    ),
    make_option(
        opt_str = c("-x", "--external"), 
        action  =  'store', 
        type  =  "character", 
        default = 'FALSE',
        help = "A space delimeted file with a single line containing external parameters from bowtie [default %default]",
        dest = "externalParameters"
    ),
    make_option(
        opt_str = c("-s", "--sufix"), 
        type = "character", 
        default = "genome",
        help = "Write Ebwt data to files with this indexDir/basename [default %default]",
        dest = "indexBaseName"
    ),
    make_option(
        opt_str = c("-S", "--samtools"), 
        action = "store_true", 
        default = FALSE,
        help = "Use this option if you want to convert the SAM files to sorted BAM. samtools is required [%default]",
        dest = "samtools"
    ),
    make_option(
        opt_str = c("-D", "--delete"), 
        action = "store_true", default = FALSE,
        help = "Use this option if you want to delete the SAM files after convert to sorted BAM. [%default]",
        dest = "deleteSAMfiles"
    ),
    # make_option(
    #     opt_str = c("-e", "--extractFolder"), 
    #     type = "character", 
    #     default = "03-UnmappedReadsHISAT2",
    #     help = "Directory to store the ummapped reads [default %default]",
    #     dest = "extractedFolder"
    # ),
    make_option(
        opt_str = c("-u", "--unmapped"), 
        action = "store_true", 
        default = FALSE,
        help = "Run samtools to extract unmapped reads from bam or sam files. [%default]",
        dest = "unmapped"
    )
    # make_option(
    #     opt_str = c("-g", "--gtfTargets"), 
    #     type = "character", 
    #     default = "gtf_targets.gtf",
    #     help = "Path to the gtf file [target gtf] to run mapping against. If would like to run without gtf file, -g option is not required [default %default]",
    #     dest = "gtfTarget"
    #     ),
    # make_option(),
    # make_option(),
    # make_option(),
    # make_option(),
    # make_option(),
    # make_option()
    
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list, description =  paste('Author: OLIVEIRA, H.C.', 'Version: 0.3.5', 'E-mail: hanielcedraz@gmail.com', sep = "\n", collapse = '\n')))



write(glue("\n\n {str_dup('=', 100)} \n\n {str_dup(' ', 40)} Mapping reads using {opt$mappingProgram} \n\n {str_dup('=', 100)} \n\n"), stdout())

write(glue("\n\n {str_dup('\n', 5)} \n\n"), stdout())
# write(glue("\n\n {str_dup('*', 40)} \n\n"), stdout())
# write(glue("\n\n {str_dup('*', 40)} \n\n"), stdout())
# write(glue("\n\n {str_dup('*', 40)} \n\n"), stdout())



if (is.null(opt$mappingFolder)) {
    if (opt$mappingProgram == "bowtie") {
        opt$mappingFolder = "02-MappedReadsBowtie"
        
    } else if (opt$mappingProgram == "bowtie2"){
        opt$mappingFolder = "02-MappedReadsBowtie2"
        
    } else if (opt$mappingProgram == "bwa"){
        opt$mappingFolder = "02-MappedReadsBWA"
        
    } else {
        stop("Error in setting mapping folder")
    }
}



if (!file.exists(file.path(opt$mappingFolder))) {
    dir.create(file.path(opt$mappingFolder), recursive = TRUE, showWarnings = FALSE)
    
}


########################################
### CHECKING IF BOWTIE EXIST IN THE SYSTEM
########################################



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
#         write(paste("umi_tools was installed in", pwdInstal, "UMI-Tools/bin/"), stderr())
#         umi_tools <- paste(pwdInstal, "UMI-Tools/bin/umi_tools")
#     }
# }


# if (detectCores() < opt$procs) {
#     write(paste("number of cores specified (", opt$procs,") is greater than the number of cores available (",detectCores(),")"), stdout())
#     paste('Using ', detectCores(), 'threads')
# }


externalPar <- opt$externalParameters
if (file.exists(externalPar)) {
    con = file(externalPar, open = "r")
    line = readLines(con, warn = FALSE, ok = TRUE)
    write(
        glue(
            "\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Using external commands: {line} \n\n {str_dup('-', 100)} \n\n"
        ), 
        stdout()
    )
}


# verify if sample_file exist
if (!file.exists(opt$samplesFile)) {
    stop(paste("Sample file", opt$samplesFile, "does not exist\n"))
}


## loadSampleFile
cat("\n\n")

write(
    glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Samples to process \n\n"), 
    stdout()
)
#opt$libraryType <- "pairEnd"

# samples <- baqcomPackage::loadSamplesFile(
#     file = "samplesSingle.txt",
#     reads_folder = "01-CleanedReadsSingle/",
#     column = opt$samplesColumn, libraryType = "singleEnd"
# )
# # 
# baqcomPackage::createSampleList(
#     samples = samples,
#     reads_folder = "01-CleanedReadsSingle/",
#     column = opt$samplesColumn,
#     libraryType = "singleEnd",
#     program = opt$mappingProgram
# )

#opt$samplesFile <- "/Users/haniel/OneDrive/posDoc/miRbaqcom/samplesReal.txt"
#opt$cleanedFolder <- "/Users/haniel/OneDrive/posDoc/miRbaqcom/01-CleanedReadsReal/"
samples <- baqcomPackage::loadSamplesFile(
    file = opt$samplesFile, 
    reads_folder = opt$cleanedFolder, 
    column = opt$samplesColumn, 
    libraryType = opt$libraryType
) 


cat("\n\n")

print(samples)
write(glue("\n\n {str_dup('-', 100)} \n\n"), stdout())
write(glue("\n\n {str_dup('\n', 3)} \n\n"), stdout())

write(
    glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 30)} Numbers of threads to be used \n\n "), 
    stdout()
)
procs <- baqcomPackage::prepareCore(nThreads = opt$procs)
#write(glue("\n\n {str_dup('\n', 5)} \n\n"), stdout())
write(glue("\n\n {str_dup('-', 100)} \n\n"), stdout())
write(glue("\n\n {str_dup('\n', 3)} \n\n"), stdout())
write(
    glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Parameters \n\n "), 
    stdout()
)

write(
    glue("\n Samples file: {opt$samplesFile} \n Input folder: {opt$cleanedFolder} \n library type: {opt$libraryType} \n Mapping program: {opt$mappingProgram} \n Output folder: {opt$mappingFolder} \n Samples to process in parallel: {opt$sampleToprocs} samples \n\n {str_dup('-', 100)}"),
    stdout()
)




# ext <- unique(file_ext(dir(file.path(opt$cleanedFolder), pattern = "gz")))
# if (length(ext) == 0) {
#     write(paste("Cannot locate fastq or sff file in folder", reads_folder, "\n"), stderr())
#     stop()
# }

if (opt$mappingProgram == "bowtie") {
    
    # write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Mapping started  using {opt$libraryType} reads \n\n {str_dup('-', 100)} \n\n"), stdout())
    
    
    #write(glue("\n\n Uncompressing files............ \n\n "), stdout())
    #if (!all(file.exists(list_files_with_exts(index_Folder, exts = ".gz"))))
    
    filetype <- function(path){
        f = file(path)
        ext = summary(f)$class
        close.connection(f)
        ext
    }
    #opt$cleanedFolder <- "01-CleanedReadsPair/"
    # file("00-FastqSingle/SRR13450790_SE_001.fastq.gz")
    allFiles <- list.files(opt$cleanedFolder, full.names = TRUE)
    # files <- lapply(allFiles, function(x) filetype(x))
    # filetype("01-CleanedReadsSingle/SRR13450790_trim_SE.fastq.gz")
    
    files <- lapply(allFiles, function(x) stringr::str_detect(x, pattern = "gz"))
    #files
    files <- unlist(files)
    #files
    #stringr::str_detect(allFiles, pattern = "gz")
    
    print(files)
    if (all(files)) {
        write(glue("\n\n Bowtie does not allow to use gz files, so it needs to be uncompressed before running bowtie \n\n Uncompressing files............ \n\n"), stdout())
        
        unpigz <- glue::glue("unpigz -p {procs} {allFiles}")
        
        write(glue("\n\n {str_dup('-', max(str_count(unpigz)))} \n \n"), stdout())
        for (i in 1:length(unpigz)) {
            print(unpigz[i])
            system(paste(unpigz[i]))
        }
        write(glue("\n\n {str_dup('-', max(str_count(unpigz)))} \n \n"), stdout())
    }
    
    
    
    #"SRR13450790_SE_001.fastq"
    
    # samples <- samples %>%
    #     mutate(Read_1 = str_remove(Read_1, ".gz"))
    # 
    # #opt$cleanedFolder <- "00-Fastq/"
    # MappingQuery <- baqcomPackage::createSampleList(
    #     samples = samples,
    #     reads_folder = opt$cleanedFolder,
    #     column = opt$samplesColumn,
    #     libraryType = opt$libraryType,program = opt$mappingProgram,
    # )
    
}



#write(glue("\n\n {str_dup('-', 100)} \n\n"), stdout())
write(glue("\n\n {str_dup('\n', 3)} \n\n"), stdout())
write(
    glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Total of samples to process \n\n "), 
    stdout()
)



MappingQuery <- baqcomPackage::createSampleList(
    samples = samples, 
    reads_folder = opt$cleanedFolder, 
    column = opt$samplesColumn,
    libraryType = opt$libraryType, 
    program = opt$mappingProgram
)
write(
    glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} "), 
    stdout()
)
write(glue("\n\n {str_dup('\n', 3)} \n\n"), stdout())


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
#                       paste0(opt$Raw_Folder, '/',index$sampleName,'*'),
#                       #paste0(opt$Raw_Folder, '/',index$R2),
#                       ' -o ',
#                       paste0(reports,'/', beforeQC),
#                       '-t', opt$sampleToprocs
#                 )
#             )
#         })
#     }, mc.cores = opt$sampleToprocs)
#     # for (i in samples[,1]) {
#     #         system2('fastqc',
#     #                paste0(opt$Raw_Folder, '/', i,'*', ' -o ', paste0(reports,'/', beforeQC), ' -t ', opt$sampleToprocs))
#     #     }
# }






## create output folder
# output_Folder <- opt$output
# if (opt$umiCommand == "extract") {
#     extractedFolder <- paste0(output_Folder, "/00-extracted")
#     if (!file.exists(file.path(extractedFolder))) {
#         dir.create(file.path(extractedFolder), recursive = TRUE, showWarnings = FALSE)
#     }
# } else if (opt$umiCommand == "dedup") {
#     dedupFolder <- paste0(output_Folder, "/01-deduplicated")
#     if (!file.exists(file.path(dedupFolder))) {
#         dir.create(file.path(dedupFolder), recursive = TRUE, showWarnings = FALSE)
#     }
# } else if (opt$umiCommand == "group") {
#     dedupFolder <- paste0(output_Folder, "/02-grouped")
#     if (!file.exists(file.path(dedupFolder))) {
#         dir.create(file.path(dedupFolder), recursive = TRUE, showWarnings = FALSE)
#     }
# }







####################
### GENOME GENERATE
####################

#gtf <- if(file.exists(opt$gtfTarget)){paste('--sjdbGTFfile', opt$gtfTarget)}
# index_Folder <- "/Users/haniel/Documents/BAQCOM/examples/genome/index_BOWTIE"
#opt$mappingTarget <- "/Users/haniel/Documents/BAQCOM/examples/genome/Sus.Scrofa.chr1.genome.dna.toplevel.fa"
#index_Folder <- paste0(dirname(opt$mappingTarget), '/', 'index_', toupper(opt$mappingProgram), '/')


index_Folder <- dirname(opt$mappingTarget)
# if (!file.exists(file.path(paste(index_Folder, '/', 'Genome', sep = ''))))
# {
#     dir.create(file.path(index_Folder), recursive = TRUE, showWarnings = FALSE)
# }

if (!file.exists(index_Folder)) {
    dir.create(file.path(index_Folder), recursive = TRUE, showWarnings = FALSE)
}

#bowtie <- "/Users/haniel/miniconda3/bin/bowtie"
#"bowtie-build" <- "/Users/haniel/miniconda3/bin/bowtie-build"
indexBuiding <- function(program, mappingTarget, index_Folder) {
    
    if (program == "bowtie") {
        try({
            system(
                paste(
                    "bowtie-build", 
                    mappingTarget,
                    paste0(index_Folder, "/", opt$indexBaseName)
                )
            )
        })
        
    } else if (program == "bowtie2") {
        try({
            system(
                paste(
                    "bowtie-build", 
                    mappingTarget,
                    index_Folder
                )
            )
        })
        
    } else if (program == "bwa") {
        try({
            system(
                glue(
                    "bwa index", 
                    mappingTarget,
                    index_Folder,
                    .sep = " "
                )
            )
        })
    } 
    
    
}

#index_genom <- star.index.function()
#tools::list_files_with_exts(index_Folder, exts = "ebwt")
# if (!all(file.exists(list_files_with_exts(index_Folder, exts = "ebwt")))) {
#     
#     index_genom <- indexBuiding(program = opt$mappingProgram, opt$mappingTarget, index_Folder)
# }

userInput <- function(question) {
    cat(question)
    con <- file("stdin")
    on.exit(close(con))
    n <- readLines(con, n = 1)
    return(n)
}


# if (!(userInput("Would you like to delete and re-run index generation? (yes or no) ") %in% c("yes", "no"))) {
#   cat('\n')
#   write(paste('May have a mistake with the argument in -s parameter. Please verify if the argument is written in the right way'), stderr())
#   stop()
# }

if (!all(file.exists(list_files_with_exts(index_Folder, exts = "ebwt")))) {
    write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Buiding genome started \n\n {str_dup('-', 100)} \n\n"), stdout())
    
    index_genom <- indexBuiding(program = opt$mappingProgram, opt$mappingTarget, index_Folder)
    
    write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Buiding genome finished \n\n {str_dup('-', 100)} \n\n"), stdout())
}

if (opt$indexBuild) {
    if (!all(file.exists(list_files_with_exts(index_Folder, exts = "ebwt")))) {
        write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Buiding genome started \n\n {str_dup('-', 100)} \n\n"), stdout())
        
        index_genom <- indexBuiding(program = opt$mappingProgram, opt$mappingTarget, index_Folder)
        
        write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Buiding genome finished \n\n {str_dup('-', 100)} \n\n"), stdout())
    } else{
        write(paste("Index genome files already exists."), stderr())
        repeat {
            inp <- userInput("Would you like to delete and re-run index generation? ([yes] or no) ")
            #imp <- "yes"
            if (inp %in% c("yes", "no", "", "Y", "y", "N", "n")) {
                write(glue("\n\n Buiding genome skiped \n\n"), stdout())
                break()
                
            } else {
                write("Specify 'yes' or 'no'", stderr())
            }
        }
        if (any(inp %in% c("yes", "", "Y", "y"))) {
            write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Buiding genome started \n\n {str_dup('-', 100)} \n\n"), stdout())
            
            index_genom <- indexBuiding(program = opt$mappingProgram, opt$mappingTarget, index_Folder)
            
            write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Buiding genome finished \n\n {str_dup('-', 100)} \n\n"), stdout())
        }
    }
    
}




#Mapping
#pigz <- system('which pigz 2> /dev/null', ignore.stdout = TRUE, ignore.stderr = TRUE)
# bowtie -S -p 8 /Users/haniel/Documents/BAQCOM/examples/genome/genome 00-Fastq/SRR13450790_SE_001.fastq > bowtieTest

#print(MappingQuery)

indexFiles <- paste0(index_Folder, "/", opt$indexBaseName)


write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 30)} Mapping started  using {opt$libraryType} reads \n\n {str_dup('-', 100)} \n\n"), stdout())

if (opt$mappingProgram == "bowtie") {
    
    # write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Mapping started  using {opt$libraryType} reads \n\n {str_dup('-', 100)} \n\n"), stdout())
    
    
    #write(glue("\n\n Uncompressing files............ \n\n "), stdout())
    #if (!all(file.exists(list_files_with_exts(index_Folder, exts = ".gz"))))
    # 
    # filetype <- function(path){
    #     f = file(path)
    #     ext = summary(f)$class
    #     close.connection(f)
    #     ext
    # }
    # #opt$cleanedFolder <- "01-CleanedReadsPair/"
    # # file("00-FastqSingle/SRR13450790_SE_001.fastq.gz")
    # allFiles <- list.files(opt$cleanedFolder)
    # # files <- lapply(allFiles, function(x) filetype(x))
    # # filetype("01-CleanedReadsSingle/SRR13450790_trim_SE.fastq.gz")
    # 
    # files <- lapply(allFiles, function(x) stringr::str_detect(x, pattern = "gz"))
    # #files
    # files <- unlist(files)
    # #files
    # #stringr::str_detect(allFiles, pattern = "gz")
    # 
    # 
    # if (all(files)) {
    #     write(glue("\n\n Bowtie does not allow to use gz files, so it needs to be uncompressed before running bowtie \n\n Uncompressing files............ \n\n"), stdout())
    #     
    #     unpigz <- glue::glue("unpigz -p {procs} {allFiles} ")
    #     
    #     write(glue("\n\n {str_dup('-', max(str_count(unpigz)))} \n \n"), stdout())
    #     for (i in 1:length(unpigz)) {
    #         print(unpigz[i])
    #         system2(unpigz[i])
    #     }
    #     write(glue("\n\n {str_dup('-', max(str_count(unpigz)))} \n \n"), stdout())
    # }
    # 
    # 
    # 
    # #"SRR13450790_SE_001.fastq"
    # 
    # samples <- samples %>%
    #     mutate(Read_1 = str_remove(Read_1, ".gz"))
    # 
    # #opt$cleanedFolder <- "00-Fastq/"
    # MappingQuery <- baqcomPackage::createSampleList(
    #     samples = samples,
    #     reads_folder = opt$cleanedFolder,
    #     column = opt$samplesColumn,
    #     libraryType = opt$libraryType,program = opt$mappingProgram,
    # )
    # 
    # 
    # 
    # #print(MappingQuery)
    # if (opt$libraryType == "singleEnd") {
    #     
        
        
        # write(glue("\n\n {str_dup('-', 100)} \n\n Starting Bowtie Mapping Single-End \n\n {str_dup('-', 100)}"), stdout())
        
        
        bowtieSingle <- mclapply(MappingQuery, function(index){
            # printing used command in terminal
            write(paste("Command:",
                        paste(
                            "bowtie",
                            "-S",
                            paste("-p", procs),
                            indexFiles,
                            paste0(index$SE),
                            if (file.exists(externalPar)) line,
                            paste0("> ", opt$mappingFolder, "/", index$sampleName, "_unsorted_sample.sam")
                        )
            ), stdout())
            
            try({
                system(
                    paste(
                        "bowtie",
                        "-S",
                        paste("-p", procs),
                        indexFiles,
                        paste0(index$SE),
                        if (file.exists(externalPar)) line,
                        paste0("> ", opt$mappingFolder, "/", index$sampleName, "_unsorted_sample.sam")
                    )
                )
            })
        }, mc.cores = opt$sampleToprocs)
        
        if (!all(sapply(bowtieSingle , "==", 0L))) {
            stop(paste("Something went wrong with Bowtie. Some jobs failed"))
        }
        
        
    } else if (opt$libraryType == "pairEnd") {
        bowtiePair <- mclapply(MappingQuery, function(index){
            # printing used command in terminal
            write(paste("Command:",
                        paste(
                            "bowtie",
                            "-S",
                            paste("-p", procs),
                            indexFiles,
                            paste0(index$PE1),
                            paste0(index$PE2),
                            if (file.exists(externalPar)) line,
                            paste0("> ", opt$mappingFolder, "/", index$sampleName, "_unsorted_sample.sam")
                        )
            ), stdout())
            try({
                system(
                    paste(
                        "bowtie",
                        "-S",
                        paste("-p", procs),
                        indexFiles,
                        paste0(index$PE1),
                        paste0(index$PE2),
                        if (file.exists(externalPar)) line,
                        paste0("> ", opt$mappingFolder, "/", index$sampleName, "_unsorted_sample.sam")
                    )
                )
            })
        }, mc.cores = opt$sampleToprocs)
        
        if (!all(sapply(bowtiePair , "==", 0L))) {
            stop(paste("Something went wrong with Bowtie. Some jobs failed"))
        }
    }
    
    write(glue("\n\n Bowtie mapping has finished \n\n "), stdout())
    
    files <- list.files(opt$cleanedFolder, full.names = TRUE)
    if (!all(lapply(files, function(x) filetype(x)) == "gzfile")) {
        write(glue("\n\n Compressing files............ \n\n "), stdout())
        system(
            paste("pigz", paste("-p", procs), paste0(opt$cleanedFolder, "/*")
            )
        )
    }
    
} else if (opt$mappingProgram == "bowtie2") {
    if (opt$libraryType == "singleEnd") {
        bowtie2Pair <- mclapply(MappingQuery, function(index){
            try({
                system(
                    paste(
                        "umi_tools",
                        "extract",
                        paste("-I", index$R1),
                        paste("--bc-pattern", "NNNXXXXNN", sep = "="),
                        paste("--read2-in", index$R2, sep = "="),
                        paste("--stdout", paste0(extractedFolder, "/", index$sampleName, "_extracted_R1.fastq.gz"), sep = "="),
                        paste("--read2-out", paste0(extractedFolder, "/", index$sampleName, "_extracted_R2.fastq.gz"), sep = "=")
                    )
                )
            })
        }, mc.cores = opt$sampleToprocs)
        
        if (!all(sapply(bowtie2Pair , "==", 0L))) {
            stop(paste("Something went wrong with Bowtie2. Some jobs failed"))
        }
        
    } else if (opt$libraryType == "pairEnd") {
        bowtie2Single <- mclapply(MappingQuery, function(index){
            try({
                system(
                    paste(
                        "umi_tools",
                        "extract",
                        paste("-I", index$R1),
                        paste("--bc-pattern", "NNNXXXXNN", sep = "="),
                        paste("--read2-in", index$R2, sep = "="),
                        paste("--stdout", paste0(extractedFolder, "/", index$sampleName, "_extracted_R1.fastq.gz"), sep = "="),
                        paste("--read2-out", paste0(extractedFolder, "/", index$sampleName, "_extracted_R2.fastq.gz"), sep = "=")
                    )
                )
            })
        }, mc.cores = opt$sampleToprocs)
        
        if (!all(sapply(bowtie2Single , "==", 0L))) {
            stop(paste("Something went wrong with Bowtie2. Some jobs failed"))
        }
    }
    
} else if (opt$mappingProgram == "bwa") {
    
    
    if (opt$libraryType == "singleEnd") {
        bwaSingle <- mclapply(MappingQuery, function(index){
            write(paste('Starting Single-End Mapping sample', index$sampleName), stderr())
            try({
                system(paste('bwa mem',
                             '-t', ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs)),
                             opt$mappingTarget,
                             #paste0(index_Folder,index_names),
                             paste0(index$SE, " > ", opt$mappingFolder, '/', index$sampleName, '_unsorted_sample.sam'),
                             if (file.exists(externalPar)) line))})
        }, mc.cores = opt$sampleToprocs
        )
        
        
        if (!all(sapply(bwaSingle, "==", 0L))) {
            stop(paste("Something went wrong with BWA mapping. Some jobs failed"))
        }
        
    } else if (opt$libraryType == "pairEnd") {
        bwaPair <- mclapply(MappingQuery, function(index){
            write(paste('Starting Single-End Mapping sample', index$sampleName), stderr())
            try({
                system(paste('bwa mem',
                             '-t', ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs)),
                             opt$mappingTarget,
                             paste0(index$PE1),
                             paste0(index$PE2),
                             paste0(">", opt$mappingFolder, '/', index$sampleName, '_unsorted_sample.sam'),
                             if (file.exists(externalPar)) line))})
        }, mc.cores = opt$sampleToprocs
        )
        
        
        if (!all(sapply(bwaPair, "==", 0L))) {
            stop(paste("Something went wrong with BWA mapping. Some jobs failed"))
        }
    }
}




santools.map <- createSampleList(
    samples = samples,
    reads_folder = opt$mappingFolder,
    column = opt$samplesColumn, fileType = "sam",
    libraryType = opt$libraryType,
    program = "samtools"
)


## run samtools
if (opt$samtools) {
    write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 20)} Sorting genes by coordenates and generating stats\n\n {str_dup('-', 100)} \n\n"), stdout())
    samtools_out <- mclapply(santools.map, function(index){
        #dir.create(file.path(opt$mappingFolder, index$sampleName))
        try({
            system(
                paste(
                    'samtools',
                    'sort',
                    paste('--threads', procs),
                    paste0(index$unsorted_sample, collapse = ","),
                    '>', paste0(opt$mappingFolder,'/', index$sampleName, '_sam_sorted_pos.bam')
                )
            )
            
            # system(
            #     paste(
            #         'samtools',
            #         'depth -a',
            #         paste('--threads', procs),
            #         paste0(index$unsorted_sample, collapse = ","),
            #         paste("| awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'"),
            #         paste0(" > ", opt$mappingFolder,'/', index$sampleName, '_sam_sorted_pos.depth')
            #     )
            # )
            
            #samtools depth -a file.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'
            
            system(
                paste(
                    'samtools',
                    'flagstat',
                    paste('--threads', procs),
                    paste0(index$unsorted_sample, collapse = ","),
                    '>', paste0(opt$mappingFolder,'/', index$sampleName, '_sam_sorted_pos.flagstat')
                )
            )
            
            
            
        })
    }, mc.cores = procs)
    
    if (!all(sapply(samtools_out, "==", FALSE))){
        write(paste("Something went wrong with samtools processing some jobs failed"),stderr())
        stop()
    }
    
    
    
    # samtools.run <- mclapply(santools.map, function(index){
    #     write(paste('Starting convert sam to bam with samtools:', index$sampleName), stderr())
    #     try({
    #         system(paste('samtools',
    #                      'sort',
    #                      paste('--threads', procs),
    #                      paste0(index$unsorted_sample, collapse = ","),
    #                      '>', paste0(opt$mappingFolder,'/', index$sampleName, '_sam_sorted_pos.bam')))})
    # }, mc.cores = opt$sampleToprocs
    # )
    # 
    # 
    # if (!all(sapply(samtools.run, "==", 0L))) {
    #     stop(paste("Something went wrong with SAMTOOLS. Some jobs failed"))
    #     
    # }
    
} else {
    
    
    write(glue("\n\n {str_dup('-', 100)} \n\n {str_dup(' ', 40)} Generating stats\n\n {str_dup('-', 100)} \n\n"), stdout())
    
    samtools_out <- mclapply(santools.map, function(index){
        #dir.create(file.path(opt$mappingFolder, index$sampleName))
        try({
            # system(
            #     paste(
            #         'samtools',
            #         'depth -a',
            #         paste('--threads', procs),
            #         paste0(index$unsorted_sample, collapse = ","),
            #         paste("| awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'"),
            #         paste0(" > ", opt$mappingFolder,'/', index$sampleName, '_sam_sorted_pos.depth')
            #     )
            # )
            
            #samtools depth -a file.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'
            
            system(
                paste(
                    'samtools',
                    'flagstat',
                    paste('--threads', procs),
                    paste0(index$unsorted_sample, collapse = ","),
                    '>', paste0(opt$mappingFolder,'/', index$sampleName, '_sam_sorted_pos.flagstat')
                )
            )
            
            
            
            
        })
    }, mc.cores = procs)
    
    if (!all(sapply(samtools_out, "==", FALSE))){
        write(paste("Something went wrong with samtools processing some jobs failed"),stderr())
        stop()
    }
}




# creating extracted_Folder
if (opt$unmapped) {
    if (!opt$samtools) {
        extracted_Folder <- opt$extractedFolder
        if (!file.exists(file.path(extracted_Folder))) dir.create(file.path(extracted_Folder), recursive = TRUE, showWarnings = FALSE)
        samtools.ummaped <- mclapply(santools.map, function(index){
            write(paste('Starting extract ummapped reads from sample', index$sampleName), stderr())
            try({
                system(paste('samtools',
                             'view',
                             '--threads', ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs)),
                             '-b',
                             '-f',
                             4,
                             paste0(index$unsorted_sample),
                             '>', paste0(opt$extractedFolder,'/', index$sampleName, '_unmapped_unsorted_pos.bam')))
                system(paste('samtools',
                             'bam2fq',
                             paste0(opt$extractedFolder,'/', index$sampleName, '_unmapped_unsorted_pos.bam'),
                             '>', paste0(opt$extractedFolder,'/', index$sampleName, '_unmapped_unsorted_pos.fastq')
                             
                ))
                unlink(paste0(opt$extractedFolder,'/', index$sampleName, '_unmapped_unsorted_pos.bam'))})
        }, mc.cores = opt$sampleToprocs
        )
    } else if (opt$samtools) {
        santools.map <- samtoolsList(samples, opt$inputFolder, opt$samplesColumn)
        extracted_Folder <- opt$extractedFolder
        if (!file.exists(file.path(extracted_Folder))) dir.create(file.path(extracted_Folder), recursive = TRUE, showWarnings = FALSE)
        #
        samtools.ummaped <- mclapply(santools.map, function(index){
            write(paste('Starting extract ummapped reads from sample', index$sampleName), stderr())
            try({
                system(paste('samtools',
                             'view',
                             '--threads', ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs)),
                             '-b',
                             '-f',
                             4,
                             paste0(opt$mappingFolder,'/',index$sampleName,'_sam_sorted_pos.bam'),
                             '>', paste0(opt$extractedFolder,'/', index$sampleName, '_unmapped_sorted_pos.nam')))
                system(paste('samtools',
                             'bam2fq',
                             paste0(opt$extractedFolder,'/', index$sampleName, '_unmapped_sorted_pos.bam'),
                             '>', paste0(opt$extractedFolder,'/', index$sampleName, '_unmapped_sorted_pos.fastq')
                             
                ))
                unlink(paste0(opt$extractedFolder,'/', index$sampleName, '_unmapped_unsorted_pos.bam'))})
        }, mc.cores = opt$sampleToprocs
        )
        
        
        if (!all(sapply(samtools.run, "==", 0L))) {
            write(paste("Something went wrong with SAMTOOLS. Some jobs failed"),stderr())
            stop()
        }
    }
}



if (opt$deleteSAMfiles) {
    unlink(dir(path = file.path(opt$mappingFolder), recursive = TRUE, pattern = ".sam$", full.names = TRUE))
}


#opt$mappingFolder <- "02-MappedReadsBWA"
# files <- list.files(opt$mappingFolder, pattern = ".flagstat", full.names = TRUE)
# getwd()
# stats <- lapply(files, function(x) readLines(x))
# print(stats)

#setwd("/Users/haniel/OneDrive/posDoc/miRbaqcom/02-MappedReadsBWA/")
#x <- readr::read_log("HE20-100K_sam_sorted_pos.flagstat", show_col_types = FALSE, trim_ws = TRUE)


TidyTable <- function(x) {
    final <- tibble::tibble(
        "in total (QC-passed reads + QC-failed reads)" = x[1,1],
        'primary' = x[2,1],
        'supplementary' = x[3,1],
        'duplicates' = x[4,1],
        'primary duplicates' = x[5,1],
        'mapped' = glue::glue("{x[7,1]} {x[7,5]})"),
        'primary mapped' = glue::glue("{x[7,1]} {x[8,6]})"),
        'paired in sequencing' = x[8,1],
        'read1' = x[9,1],
        'read2' = x[10,1],
        'properly paired' = glue::glue("{x[11,1]} {x[8,6]})"),
        'with itself and mate mapped' = x[12,1],
        'singletons' = glue::glue("{x[13,1]} {x[8,6]})")
    )
    return(final)
}

report_sample <- list()
#opt$mappingFolder <- "02-MappedReadsBWA/"
#print(getwd())
statFiles <- list.files(opt$mappingFolder, pattern = ".flagstat", full.names = TRUE)
for (i in statFiles) { # change this to your "samples"
    #print(i)
    suppressWarnings(report_sample[[i]] <- readr::read_log(
        i, 
        show_col_types = FALSE, 
        col_names = c(
            "in total (QC-passed reads + QC-failed reads)",
            'primary',
            'supplementary',
            'duplicates',
            'primary duplicates',
            'mapped',
            'primary mapped',
            'paired in sequencing',
            'read1',
            'read2',
            'properly paired',
            'with itself and mate mapped',
            'singletons'
        )
    )
    )
}

#list.files("/Users/haniel/OneDrive/posDoc/miRbaqcom/02-MappedReadsBWA/", pattern = ".flagstat")
df <- lapply(report_sample, FUN = function(x) TidyTable(x))
final_df <- do.call("rbind", df)
report <- "report"
if (!file.exists(report)) {
    dir.create(report, recursive = TRUE)
}
write.table(final_df, paste0(report,"/report_mapping", opt$mappingProgram, ".log"))

print(data.table::as.data.table(final_df))




# # 
# #####################################################
# sapply(samples, function(tgt) {
#     print(paste("Generating output for target:", tgt[1]))})
# files <- list.files(opt$mappingFolder, pattern = ".sam")
# 
# file.path()
# samplesID <- samples[,1]
# # write out index tables
# targetTables <- sapply(samples[,opt$samplesColumn], function(tgt) {
#     print(paste("Generating output for target:",tgt[1]))
#     filesToRead <- unlist(file.path(
#         opt$mappingFolder,
#         paste0(tgt[1], "_sam_sorted_pos.flagstat")
#     )
#     
#     
#     )
#     
#     
#     #	filesToRead <- unlist(sapply(file.path(opt$mappingFolder,unique(samples[,opt$samplesColumn])),dir,pattern=paste(tgt[1],"idxstats",sep="."),full.names=TRUE))
#     info <- read.table(filesToRead[1])[,1:2]
#     colnames(info) <- c("SequenceID","SequenceLength")
#     data <- sapply(filesToRead,function(file){
#         tb <- read.table(file)
#         values <- rowSums(as.matrix(tb[,3:4]))
#         values
#     })
#     colnames(data) <- basename(colnames(data))
#     freq <- round(sweep(data,2,colSums(data),"/")*100,3)
#     write.table(cbind(info,data),file.path(opt$mappingFolder,paste(tgt[1],"summary","reads","txt",sep=".")),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
#     write.table(cbind(info,freq),file.path(opt$mappingFolder,paste(tgt[1],"summary","proportions","txt",sep=".")),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
#     pmapped <- 100-tail(freq,1)
#     #	as.vector(pmapped)
#     names(pmapped) <- colnames(freq)
#     round(pmapped,3)
# })
# 
# targetTables <- data.frame(targetTables)
# colnames(targetTables) <- sapply(targets,"[[",1L)
# 
# ### simple assign by most on target
# targetTables <- data.frame(ID=rownames(targetTables),targetTables,assign=colnames(targetTables)[apply(targetTables,1,which.max)])
# write.table(targetTables,file.path(opt$mappingFolder,"SummarySample2Targets.txt"),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
# 
# cat('\n')
# 
# # Creating samples report
# if (!opt$singleEnd) {
#     TidyTable <- function(x) {
#         final <- data.frame('Input_Read_Pairs' = x[1,2],
#                                          header = F, as.is = T, fill = TRUE, sep = ':', text = TRUE)
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
# cat('\n')
# write(paste('How to cite:', sep = '\n', collapse = '\n', "Please, visit https://github.com/hanielcedraz/BAQCOM/blob/master/how_to_cite.txt", "or see the file 'how_to_cite.txt'"), stderr())
# 

# system(
#     glue("conda deactivate")
# )

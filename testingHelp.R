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
        dest = "bcPattern", callback = "function(--extractMethod regex){} "
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
        dest = "output"
        #metavar = "output_folder"
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
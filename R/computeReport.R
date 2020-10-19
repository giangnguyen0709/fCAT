#' This function arrange the completeness output of the interested 
#' genome into the corresponding folder in the output folder
#' 
#' @param completenessOuput a list contains the completeness output. The full
#' table, missing table and the ignored table
#' @param output The directory which contains the output directory
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param genomeName the name of the interested genome
#' @param scoreMode the mode determines the method to scoring the founded
#' ortholog and how to classify them. Choices: 1, 2, 3, "len"
#' @export

arrangeCompletenessOutput <- function(
    completenessOutput, output, coreSet, genomeName, scoreMode
) {
    if (!endsWith(output, "/")) {
        output <- paste(output, "/", sep = "")
    }
    outDir <- paste(
        output, "fcatOutput", "/", coreSet, "/", genomeName, sep = ""
    )
    if (scoreMode == "len") {
        subScoreModeDir <- paste(outDir, "/", "mode_len", sep = "")
    } else {
        subScoreModeDir <- paste(
            outDir, "/", "mode", as.character(scoreMode), sep = "")
    }

    if (!dir.exists(subScoreModeDir)) {
        dir.create(subScoreModeDir, recursive = TRUE)
    }
    
    ## Printing summary text
    summaryTable <- translateReport(
        genomeName, completenessOutput[[1]], scoreMode
    )
    summaryFile <- paste(subScoreModeDir, "/", "summary.txt", sep = "")
    
    write.table(
        summaryTable,
        summaryFile,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    
    write.table(
        completenessOutput[[1]],
        paste(subScoreModeDir, "/", "full_table.txt", sep = ""),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    
    write.table(
        completenessOutput[[2]],
        paste(outDir, "/", "missing.txt", sep = ""),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    
    write.table(
        completenessOutput[[3]],
        paste(outDir, "/", "ignored.txt", sep = ""),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}

#' The function take a path to a genome fasta file and its FAS annotation
#' file (if the annotation file is not provided, the funtion will compute it) 
#' and compute the detailed report of the completeness of the genome
#' 
#' @usage 
#' computeReport(
#'     genome, fasAnno, root, coreSet, extend = FALSE,
#'     scoreMode, priorityList = NULL, cpu, computeOri = FALSE,
#'     blastDir = NULL, weightDir = NULL, cleanup = FALSE,
#'     reFdog = FALSE, fdogDir = NULL, ppDir = NULL
#' )
#'
#' @param genome The path to the genome fasta file
#' @param fasAnno The path to the fas annotation file. It can equal NULL
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param extend The output of the function is a phylogenetic profile of the
#' interested genome. It contains 4 files, .phyloprofile, .extended.fa,
#' _reverse.domains and _forward.domains. If extend = TRUE, the files will be
#' appended into the old files in the folder output of the core directory or in
#' the inputed folder by the user with the argument ppDir. If there is no old
#' files in the folder, the output files of the function will be writen in the
#' new files.
#' @param scoreMode the mode determines the method to scoring the founded
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @param priorityList A list contains one or many genome ID of the genomes,
#' which were used to build the core set. The genome ID of this list will be
#' stored with an priority order, the tool look at into the fasta file of each
#' core group and determine with the priority order to determine the references
#' species for each core group.
#' @param cpu determines the cores that fDOG and fdogFAS will uses to be run
#' parallel
#' @param blastDir The user can replace the blast_dir folder in the core
#' directory by specifying it in this argument. By default is NULL
#' @param weightDir The user can replace the weight_dir folder in the core
#' directory by specifying it in this argument. By default is NULL
#' @param outDir The user can specify the directory to save the output report
#' file of the completeness of the interested genome by specifying the path
#' to the folder in this argument. By default is NULL
#' @param cleanup The fDOG's output is a set of phylogenetic profile of each
#' core group to the interested genome. The phylogenetic profile will be stored
#' into a folder in the core set. The function will merge all the small
#' phylogenetic profile, calculate the FAS score or length to have the whole
#' phylogenetic profile of the interested genome to the core set. This fDOG's
#' output can be reused for all score modes. When cleanup is set to TRUE, the
#' fDOG's output will not be stored to be reused but to be removed
#' @param output The directory which contains the output directory
#'
#' @return a report of the completeness of the interested genome with detailed
#' information of every core genes in the core set. Which core gene is "similar"
#' , which is "dissimilar", "duplicated" and "missing"
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' genome <- system.file("extdata", "HUMAN@9606@3.fa", package = "fCAT")
#' fasAnno <- system.file("extdata", "HUMAN@9606@3.json", package = "fCAT")
#' report <- computeReport(genome, fasAnno, coreFolder, "test",
#'     scoreMode = 2, priorityList = c("HUMAN@9606@1"), cpu = 4,
#'     output = getwd()
#' )
#' print.data.frame(report)
#' @export
computeReport <- function(
    genome, fasAnno, root, coreSet, extend = FALSE,
    scoreMode, priorityList = NULL, cpu, computeOri = FALSE,
    blastDir = NULL, weightDir = NULL, output
) {
    startTime <- Sys.time()
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    if (!endsWith(output, "/")) {
        output <- paste(output, "/", sep = "")
    }
    outDir <- paste(output, "fcatOutput", "/", coreSet, "/", sep = "")

    splited <- strsplit(genome, "/", fixed = TRUE)[[1]]
    splited <- splited[length(splited)]
    genomeName <- strsplit(splited, ".", fixed = TRUE)[[1]][1]
    
    if (scoreMode == 1) {
        ppPath <- paste(
            outDir, genomeName, "/", "phyloprofileOutput", "/", "mode1", 
            ".phyloprofile", sep = "")
    } else {
        if (scoreMode == "len") {
            ppPath <- paste(
                outDir, genomeName, "/", "phyloprofileOutput", "/", "len", 
                ".phyloprofile", sep = "")
        } else {
            ppPath <- paste(
                outDir, genomeName, "/", "phyloprofileOutput", "/", "other", 
                ".phyloprofile", sep = "")
        }
    }
    
    if (file.exists(ppPath)) {
        placeSeedCheck <- 0
        pp <- read.table(ppPath, sep = "\t", header = TRUE)
    } else {
        placeSeedCheck <- 1
        placeSeed(genome, fasAnno, root, computeOri, weightDir)
        if (scoreMode == "len") {
            pp <- runFdogBusco(
                root, coreSet, extend, scoreMode, priorityList, cpu,
                blastDir, weightDir, output
            )
        } else {
            pp <- runFdog(
                root, coreSet, extend, scoreMode, priorityList, cpu,
                blastDir, weightDir, output
            )
        }
    }
    completenessOutput <- reportSingle(
        pp, root, coreSet, scoreMode, priorityList
    )
    
    arrangeCompletenessOutput(
        completenessOutput, output, coreSet, genomeName, scoreMode
    )
    
    if (placeSeedCheck == 1) {
        unlink(
            paste(root, "query_taxon", "/", genomeName, sep = ""),
            recursive = TRUE
        )
        
        if (computeOri == FALSE) {
            if (!is.null(weightDir)) {
                if (!endsWith(weightDir, "/")) {
                    weightDir <- paste(weightDir, "/", sep = "")
                }
                file.remove(
                    paste(weightDir, genomeName, ".json", sep = "")
                )
            } else {
                file.remove(
                    paste(
                        root, "weight_dir", "/", genomeName, ".json", sep = ""
                    )
                )
            }
        }
    }
    endTime <- Sys.time()
    print("Running time:")
    print(endTime - startTime)
}

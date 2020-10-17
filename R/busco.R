#' In the score mode "busco" the tool must calculate the length of the 
#' founded ortholog sequence to classify it into "complete" or "fragmented". 
#' The function will take the phylogenetic profile and the extended fasta file 
#' of the interested genome as the input. Based on the common ID of the 
#' sequence in each file, the function will calculate the length of the 
#' sequence in the fasta file and update it into the phylogenetic profile as 
#' a new column.
#'
#' @param pp the phylogenetic profile of the interested genome in data.frame
#' @param exFasta the extended fasta file in form  of a vector, with each
#' element is a line of the file
#'
#' @return the phylogenetic profile of the interested genome with updated
#' length for each ortholog (or for each line)
#' @examples
#' ## Create pseudo data
#' geneID <- c("530670", "530730")
#' ncbiID <- c("ncbi9606", "ncbi9606")
#' orthoID <- c("530670|HUMAN@9606@3|Q16526|1", "530730|HUMAN@9606@3|P05091|1")
#' pp <- data.frame(geneID, ncbiID, orthoID)
#'
#' fasta <- c(
#'     ">530670|HUMAN@9606@3|Q16526|1",
#'     "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
#'     ">530730|HUMAN@9606@3|P05091|1",
#'     "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
#' )
#'
#' ## Updating length
#' updatedPP <- updateLength(pp, fasta)
#' print.data.frame(updatedPP)
#' @export
updateLength <- function(pp, exFasta) {
    i <- (1:nrow(pp))

    orthoLength <- unlist(
        lapply(
            i,
            function(i, exFasta, pp) {
                orthoID <- pp$orthoID[i]
                for (j in 1:length(exFasta)) {
                    if (exFasta[j] == paste(">", orthoID, sep = "")) {
                        return(nchar(exFasta[j + 1]))
                    }
                }
            },
            exFasta = exFasta,
            pp = pp
        )
    )
    pp <- cbind(pp, length = orthoLength)
    return(pp)
}

#' This function is used as a modul in the function computeReport() when 
#' using socre mode "busco". The function computeReport() take the path to the 
#' genome fasta file and the annotation file (if the user did not input it, it 
#' will be computed) of the interested genome as one of its inputs and arrange 
#' the files into the equivalent folder of the core directory. The genome fasta
#' file will be copied and saved into the query_taxon folder. A symbolic link
#' of the annotation file will be created and saved in the folder weight_dir.
#' After all the files of the interested genome were arranged into the folder,
#' the function runFdogBusco or the function runFdog will be called, depent on
#' the using score mode, the function will take the path to the core directory
#' and specify the path of the folder core_orthologs, weight_dir, query_taxon,
#' blast_dir as the input hmmpath, weightpath, searchpath, blastpath for fDOG
#' and run fDOG on this inputs to search ortholog on the interested genome for
#' the equivalent inputed core set. When the search with fDOG is finished, the
#' function will calculate the FAS score with fdogFAS or the length of each
#' ortholog sequence based on the using score mode. The function with return
#' the phylogenetic profile of the interested genome to the core set within the
#' FAS scores or the lengths
#' 
#' @usage 
#' runFdogBusco(
#'     root, coreSet, extend = FALSE, scoreMode, priorityList = NULL, cpu,
#'     blastDir = NULL, weightDir = NULL, cleanup = FALSE,
#'     reFdog = FALSE, fdogDir = NULL, ppDir = NULL)
#'
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
#' new files..
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
#' directory by specifying it in this argument
#' @param weightDir The user can replace the weight_dir folder in the core
#' directory by specifying it in this argument
#' @param cleanup The fDOG's output is a set of phylogenetic profile of each
#' core group to the interested genome. The phylogenetic profile will be stored
#' into a folder in the core set. The function will merge all the small
#' phylogenetic profile, calculate the FAS score or length to have the whole
#' phylogenetic profile of the interested genome to the core set. This fDOG's
#' output can be reused for all score modes. When cleanup is set to TRUE, the
#' fDOG's output will not be stored to be reused but to be removed
#' @param reFdog If it already exist a fDOG's output for a specific core group
#' the tool will skip this core group and go to the next core group. If reFdog
#' is set to TRUE, the tool will remove all the existed fDOG's output and rerun
#' fDOG for all core groups of the set
#' @param fdogDir Normally the fDOG's output will be stored in the folder
#' fdogout in the core directory, but the user can specify the folder for
#' fDOG's output by specify the path to it in this argument. Notice here, is
#' the fDOG's output folder will contains the subfolder, equivalent to the name
#' of the interested genome, for example, the folder can contain "HUMAN@9606@3"
#' and "AMPQU@400682@2", for a completeness checking on an interested genome,
#' which has a subfolder in the fDOG's output folder with the same name, the
#' function will look into the subfolder to find the existed fDOG's output
#' @param ppDir The user can replace the default folder output in the core
#' directory, where the phylogenetic profiles are stored by his folder. The user
#' can specify the path to his folder in this argument
#'
#' @return phylogenetic profile of the genome in data.frame
#' @examples
#' ## Take the demo data
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' genome <- system.file("extdata", "HUMAN@9606@3.fa", package = "fCAT")
#' fasAnno <- system.file("extdata", "HUMAN@9606@3.json", package = "fCAT")
#'
#' ## Place seed
#' placeSeed(genome, fasAnno, coreFolder)
#'
#' ## Test runFdogBusco
#' pp <- runFdogBusco(coreFolder, "test",
#'     extend = FALSE, scoreMode = "busco",
#'     priorityList = c("HUMAN@9606@1"), cpu = 4
#' )
#'
#' print.data.frame(pp)
#'
#' ## Delete seed
#' fastaFolder <- paste(coreFolder, "/query_taxon/HUMAN@9606@3", sep = "")
#' annoPath <- paste(coreFolder, "/weight_dir/HUMAN@9606@3.json", sep = "")
#' unlink(fastaFolder, recursive = TRUE)
#' file.remove(annoPath)
#' @export
runFdogBusco <- function(
    root, coreSet, extend = FALSE, scoreMode, priorityList = NULL, cpu,
    blastDir = NULL, weightDir = NULL, redoFdog = FALSE, output) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    if (!endsWith(output, "/")) {
        output <- paste(output, "/", sep = "")
    }
    outDir <- paste(output, "fcat_output", "/", coreSet, "/", sep = "")

    setName <- coreSet

    genomeName <- list.dirs(paste(root, "query_taxon", sep = ""),
        recursive = FALSE,
        full.names = FALSE
    )[[1]]

    hmmPath <- paste(root, "core_orthologs", "/", coreSet, sep = "")
    if (!is.null(blastDir)) {
        blastPath <- blastDir
    } else {
        blastPath <- paste(root, "blast_dir", sep = "")
    }
    searchPath <- paste(root, "query_taxon", sep = "")
    if (!is.null(weightDir)) {
        weightPath <- weightDir
    } else {
        weightPath <- paste(root, "weight_dir", sep = "")
    }
    
    outPath <- paste(outDir, "fdog_output", "/", genomeName, sep = "")
    if (!dir.exists(outPath)) {
        dir.create(outPath, recursive = TRUE)
    }
    ### - check Data - ###
    command <- paste(
        "checkData1s",
        "-g", searchPath,
        "-b", blastPath,
        "-w", weightPath
    )
    system(command)
    if (!dir.exists(outPath)) {
        dir.create(outPath, recursive = TRUE)
    }

    if (redoFdog == TRUE) {
        for (
            folder in list.dirs(
                outPath,
                recursive = TRUE, full.names = TRUE
            )
        ) {
            unlink(folder, recursive = TRUE)
        }
    }

    outputSet <- lapply(
        list.dirs(
            paste(root, "core_orthologs", "/", coreSet, sep = ""),
            recursive = FALSE,
            full.names = FALSE
        ),
        function(
            coreGene, hmmPath, blastPath, searchPath, outPath, weightPath, root,
            coreSet, scoreMode, priorityList, cpu, genomeName
        ) {
            refSpec <- getSpec(
                paste(hmmPath, "/", coreGene, "/", coreGene, ".fa", sep = ""),
                priorityList
            )
            if (is.null(refSpec)) {
                return(NULL)
            }
            if (file.exists(
                paste(
                    outPath, "/", refSpec, "/", coreGene, "/",
                    "hamstrsearch.log",
                    sep = ""
                )
            )) {
                if (
                    file.exists(
                        paste(
                            outPath, "/", refSpec, "/", coreGene, "/", coreGene,
                            ".phyloprofile",
                            sep = ""
                        )
                    )
                ) {
                    pp <- read.table(
                        paste(
                            outPath, "/", refSpec, "/", coreGene, "/", coreGene,
                            ".phyloprofile",
                            sep = ""
                        ),
                        sep = "\t",
                        header = TRUE
                    )
                    exFasta <- readLines(
                        paste(
                            outPath, "/", refSpec, "/", coreGene, "/", coreGene,
                            ".extended.fa",
                            sep = ""
                        )
                    )
                    return(list(pp, exFasta))
                } else {
                    return(NULL)
                }
            }

            command <- paste(
                "fdog.run",
                "--seqFile", paste(
                    root, "core_orthologs", "/", coreSet, "/", coreGene,
                    "/", coreGene, ".fa",
                    sep = ""
                ),
                "--seqName", coreGene,
                "--refspec", refSpec,
                "--hmmpath", hmmPath,
                "--outpath", paste(outPath, "/", refSpec, sep = ""),
                "--blastpath", blastPath,
                "--weightpath", weightPath,
                "--searchpath", searchPath,
                "--cleanup",
                "--cpu", cpu,
                "--reuseCore",
                "--checkCoorthologsRef",
                "--countercheck",
                "--fasoff"
            )
            system(command)

            if (file.exists(paste(coreGene, ".fa", sep = ""))) {
                file.remove(paste(coreGene, ".fa", sep = ""))
            }

            if (
                !file.exists(
                    paste(
                        outPath, "/", refSpec, "/", coreGene, "/", coreGene, 
                        ".phyloprofile", sep = ""
                    )
                )
            ) {
                return(NULL)
            } else {
                exFasta <- readLines(
                    paste(
                        outPath, "/", refSpec, "/", coreGene, "/", coreGene,
                        ".extended.fa",
                        sep = ""
                    )
                )
                pp <- read.table(
                    paste(
                        outPath, "/", refSpec, "/", coreGene, "/", 
                        coreGene, ".phyloprofile", sep = ""
                    ),
                    sep = "\t",
                    header = TRUE
                )

                return(list(pp, exFasta))
            }
        },
        hmmPath = hmmPath,
        blastPath = blastPath,
        searchPath = searchPath,
        weightPath = weightPath,
        outPath = outPath,
        root = root,
        coreSet = coreSet,
        priorityList = priorityList,
        cpu = cpu,
        genomeName = genomeName,
        scoreMode = scoreMode
    )

    pp <- NULL
    exFasta <- NULL
    for (output in outputSet) {
        if (is.null(pp)) {
            pp <- output[[1]]
        } else {
            pp <- rbind(pp, output[[1]])
        }
        if (is.null(exFasta)) {
            exFasta <- output[[2]]
        } else {
            exFasta <- c(exFasta, output[[2]])
        }
    }

    pp <- extractPP(pp, genomeName)
    exFasta <- extractFasta(exFasta, genomeName)

    pp <- updateLength(pp, exFasta)

    ppDir <- paste(outDir, "phyloprofile_output", "/", sep = "")

    if (!dir.exists(ppDir)) {
        dir.create(ppDir, recursive = TRUE)
    }
    
    ## Printing
    subPPDir <- paste(
        ppDir, "/", "mode_len", "/", genomeName, sep = ""
    )
    if (!dir.exists(subPPDir)) {
        dir.create(subPPDir, recursive = TRUE)
    }
    write.table(
        pp, 
        paste(subPPDir, "/", genomeName, ".phyloprofile", sep = ""),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    writeLines(
        exFasta, 
        paste(subPPDir, "/", genomeName, ".fa", sep = "")
    )

    if (extend == TRUE) {
        genomeExist <- 0
        if (
            file.exists(
                paste(
                    ppDir, folder, "/", setName, "_len.phyloprofile", sep = ""
                )
            )
        ) {
            oriPP <- read.table(
                paste(
                    ppDir, folder, "/", setName, "_len.phyloprofile", sep = ""
                ),
                sep = "\t",
                header = TRUE
            )
            genomeID <- unlist(lapply(
                oriPP$orthoID,
                function(orthoID) {
                    return(strsplit(orthoID, "|", fixed = TRUE)[[1]][2])
                }
            ))
            if (genomeName %in% genomeID) {
                genomeExist <- 1
            }
            if (genomeExist == 1) {
                while(TRUE) {
                    removeCheck <- readline(
                        prompt = "It already exists a phylogenetic profile of
                        the genome in the original phylogenetic profile.
                        Do you want to remove the old one to 
                        append the new? [y/n]:"
                    )
                    if (removeCheck == "y" || removeCheck == "n") {
                        break
                    }
                }
                if (removeCheck == "y") {
                    correctFiles(
                        paste(
                            ppDir, folder, sep = ""
                        ),
                        genomeName
                    )
                    genomeExist <- 0
                }
            }
        }
        
        if (genomeExist == 0) {
            folder <- "mode_len"
            if (
                file.exists(
                    paste(
                        ppDir, folder, "/", setName, "_len.phyloprofile", 
                        sep = ""
                    )
                )
            ) {
                oriPP <- read.table(
                    paste(
                        ppDir, folder, "/", setName, "_len.phyloprofile", 
                        sep = ""
                    ),
                    sep = "\t",
                    header = TRUE
                )
                oriPP <- rbind(oriPP, pp)
            } else {
                oriPP <- pp
            }
            
            if (
                file.exists(
                    paste(
                        ppDir, folder, "/", setName, "_len.extended.fa", 
                        sep = ""
                    )
                )
            ) {
                oriFasta <- readLines(
                    paste(
                        ppDir, folder, "/", setName, "_len.extended.fa", 
                        sep = ""
                    )
                )
                oriFasta <- c(oriFasta, exFasta)
            } else {
                oriFasta <- exFasta
            }
            
            write.table(
                oriPP,
                paste(
                    ppDir, folder, "/", setName, "_len.phyloprofile", 
                    sep = ""
                ),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
            )
            writeLines(
                oriFasta,
                paste(
                    ppDir, folder, "/´", setName, "_len.extended.fa", 
                    sep = ""
                )
            )
        }
    }
    
    ## Printing priority list
    
    table <- data.frame(
        genomeID = c(genomeName),
        ref_spec_list = c(paste(priorityList, collapse = ","))
    )
    
    priorityFile <- paste(outDir, setName, "_refspec_list", sep = "")
    
    if (!file.exists(priorityFile)) {
        write.table(
            table,
            priorityFile,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    } else {
        priorityTable <- read.table(priorityFile, header = TRUE, sep = "\t")
        if (!(genomeName %in% priorityTable$genomeID)) {
            write.table(
                table,
                priorityFile,
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE,
                quote = FALSE,
                sep = "\t"
            )
        }
    }
    
    return(pp)
}

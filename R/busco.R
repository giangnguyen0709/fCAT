#' In the scoreMode "busco" the tool must calculate the length of the founded
#' ortholog to perfome the assessment. After HaMStR founded an ortholog the
#' function will be called to update the length of the ortholog in phylogenetic
#' profile
#'
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param coreGene the ID of the core gene
#'
#' @return none
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
    ))
    pp <- cbind(pp, length = orthoLength)
    return(pp)
}

#' This function takes a path to a core set and run HaMStR to search ortholog on
#' the genome in the folder genome_dir of core set
#'
#' @param root Path to the root folder
#' @param coreSet The core set name
#' @param extend if extend=TRUE the phylogenetic profile of the genome will be
#' appended to the original phylogenetic profile
#' @param scoreMode the mode determines the way to assess the founded ortholog
#' @param priorityList the list determines the references species
#' @param cpu determines the cores that HaMStR will use
#' @param blastDir point to the user's blast_dir folder
#' @param weightDir point to the user's weight_dir folder
#' @param outDir point to the user's output folder
#' @param cleanup a logical value to decide if the fDOG 's output must be 
#' removed
#'
#' @return phylogenetic profile of the genome
#' @export
runFdogBusco <- function(
    root, coreSet, extend = FALSE, scoreMode, priorityList = NULL, cpu,
    blastDir = NULL, weightDir = NULL, outDir = NULL, cleanup = FALSE
) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    
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
    outPath <- paste(
        root, "hamstrout", "/", genomeName,
        sep = ""
    )
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
                    outPath, "/", coreGene, "/", "hamstrsearch.log", 
                    sep = ""
                )
            )) {
                if (
                    file.exists(
                        paste(
                            outPath, "/", coreGene, "/", coreGene, 
                            ".phyloprofile", sep = ""
                        )
                    )
                ) {
                    pp <- read.table(
                        paste(
                            outPath, "/", coreGene, "/", coreGene, 
                            ".phyloprofile", sep = ""
                        ),
                        sep = "\t",
                        header = TRUE
                    )
                    exFasta <- readLines(
                        paste(
                            outPath, "/", coreGene, "/", coreGene, 
                            ".extended.fa", sep = ""
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
                    "/", coreGene, ".fa", sep = ""
                ),
                "--seqName", coreGene,
                "--refspec", refSpec,
                "--hmmpath", hmmPath,
                "--outpath", outPath,
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
                        outPath, "/", coreGene, "/", coreGene, ".phyloprofile", 
                        sep = ""
                    )
                )
            ){
                return(NULL)
            } else {
                exFasta <- readLines(
                    paste(
                        outPath, "/", coreGene, "/", coreGene, 
                        ".extended.fa", sep = ""
                    )
                )
                pp <- read.table(
                    paste(
                        outPath, "/", coreGene, "/", coreGene, ".phyloprofile",
                        sep = ""
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

    if (!is.null(outDir)) {
        if (!endsWith(outDir, "/")) {
            outDir <- paste(outDir, "/", sep = "");
        }
    } else {
        outDir <- paste(
            root, "output", "/",  setName, "/", as.character(scoreMode), "/", 
            sep = ""
        )
    }
    if (!dir.exists(outDir)) {
        dir.create(outDir, recursive = TRUE)
    }
    
    if (cleanup == TRUE) {
        unlink(outPath, recursive = TRUE)
    }
    
    if (extend == TRUE) {
        if (
            file.exists(
                paste(outDir, setName, ".phyloprofile", sep = "")
            )
        ) {
            oriPP <- read.table(
                paste(outDir, setName, ".phyloprofile", sep = ""),
                sep = "\t",
                header = TRUE
            )
            oriPP <- rbind(oriPP, pp)
        } else {
            oriPP <- pp
        }
        
        if (
            file.exists(
                paste(outDir, setName, ".extended.fa", sep = "")
            )
        ) {
            oriFasta <- readLines(
                paste(outDir, setName, ".extended.fa", sep = "")
            )
            oriFasta <- c(oriFasta, exFasta)
        } else {
            oriFasta <- exFasta
        }
        
        write.table(
            oriPP,
            paste(outDir, setName, ".phyloprofile", sep = ""),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
        )
        writeLines(
            oriFasta,
            paste(outDir, setName, ".extended.fa", sep = "")
        )
    }
    return(pp)
}
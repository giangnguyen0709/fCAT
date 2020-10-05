updateExFasta <- function(
    root, coreSet,fasta, coreGene, genomeName, refSpec, scoreMode
) {
    queryFasta <- extractFasta(fasta, genomeName)
    coreFasta <- readLines(
        paste(
            root, "core_orthologs", "/", coreSet, "/", coreGene, "/", 
            coreGene, ".fa", sep = ""
        )
    )
    for (i in 1:length(coreFasta)) {
        if (i %% 2 == 1) {
            coreFasta[i] <- paste(coreFasta[i], "|1", sep = "")
        }
    }
    
    if (scoreMode != 1) {
        coreFasta <- extractFasta(coreFasta, refSpec)
    }
    
    if ((length(queryFasta) / 2) == 1) {
        queryFasta <- c(queryFasta, coreFasta)
        returnList <- list(queryFasta)
    } else {
        sequenceNumber <- length(queryFasta) / 2
        i <- 1:sequenceNumber
        returnList <- lapply(
            i,
            function(i, queryFasta, coreFasta) {
                singleFasta <- c(queryFasta[(2 * i) - 1], queryFasta[2 * i])
                singleFasta <- c(singleFasta, coreFasta)
                return(singleFasta)
            },
            queryFasta = queryFasta,
            coreFasta = coreFasta
        )
    }
    return(returnList)
}

recalculateScore <- function(pp, genomeName) {
    queryPP <- extractPP(pp, genomeName)
    genomeID <- unlist(lapply(
        pp$orthoID,
        function(orthoID) {
            return(strsplit(orthoID, "|", fixed = TRUE)[[1]][2])
        }
    ))
    singlePP <- cbind(pp, genomeID)
    singlePP <- subset(singlePP, genomeID != genomeName)
    coreGeneList <- unique(singlePP$geneID)
    fasF <- list()
    fasB <- list()
    for (gene in coreGeneList) {
        fasF[[as.character(gene)]] <- 0
        fasB[[as.character(gene)]] <- 0
        subPP <- subset(singlePP, geneID == gene)
        for (i in 1:nrow(subPP)) {
            fasF[[as.character(gene)]] <- 
                fasF[[as.character(gene)]] + subPP$FAS_F[i]
            fasB[[as.character(gene)]] <- 
                fasB[[as.character(gene)]] + subPP$FAS_B[i]
        }
        fasF[[as.character(gene)]] <- fasF[[as.character(gene)]] / nrow(subPP)
        fasB[[as.character(gene)]] <- fasB[[as.character(gene)]] / nrow(subPP)
    }
    
    for (i in 1:nrow(queryPP)) {
        queryPP[i, 4] <- fasF[[as.character(queryPP[i, 1])]]
        queryPP[i, 5] <- fasB[[as.character(queryPP[i, 1])]]
    }
    return(queryPP)
}

#' After HaMStR searched ortholog for the interested genome, it will store the
#' output in a temporary folder in the core set. The function will concanate all
#' the output files in the folder in a single files
#'
#' @param directory path to the folder, that contains the output files
#' @param genomeName the genomeID of the genome, that need to be extracted
#' @return none
#' @export
concanateFiles <- function(directory, genomeName) {
    if (!endsWith(directory, "/")) {
        directory <- paste(directory, "/", sep = "")
    }
    
    exFasta <- NULL
    domain0 <- NULL
    domain1 <- NULL
    pp <- NULL
    
    for (file in list.files(directory, full.names = TRUE, recursive = TRUE)) {
        if (endsWith(file, ".extended.fa")) {
            if (is.null(exFasta)) {
                exFasta <- readLines(file)
            } else {
                exFasta <- c(exFasta, readLines(file))
            }
        }
        
        if (endsWith(file, "_reverse.domains")) {
            domain <- try(
                read.table(file, sep = "\t", comment.char = ""), silent = TRUE)
            if (!inherits(domain, "try-error")) {
                if (is.null(domain0)) {
                    domain0 <- domain
                } else {
                    domain0 <- rbind(domain0, domain)
                }
            }
        }
        
        if (endsWith(file, "_forward.domains")) {
            domain <- try(
                read.table(file, sep = "\t", comment.char = ""), silent = TRUE)
            if (!inherits(domain, "try-error")) {
                if (is.null(domain1)) {
                    domain1 <- read.table(file, sep = "\t", comment.char = "")
                } else {
                    domain1 <- rbind(
                        domain1, read.table(
                            file, sep = "\t", comment.char = ""))
                }
            }
        }
        
        if (endsWith(file, ".phyloprofile")) {
            singlePP <- read.table(file, sep = "\t", header = TRUE)
            singlePP <- recalculateScore(singlePP, genomeName)
            if (is.null(pp)) {
                pp <- singlePP
            } else {
                pp <- rbind(pp, singlePP)
            }
        }
    }
    pp <- extractPP(pp, genomeName)
    exFasta <- extractFasta(exFasta, genomeName)
    return(list(pp, exFasta, domain0, domain1))
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
runFdog <- function(
    root, coreSet, extend = FALSE, scoreMode, priorityList = NULL, cpu,
    blastDir = NULL, weightDir = NULL, cleanup = FALSE, reFdog, fdogDir, ppDir
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
    if (!is.null(fdogDir)) {
        if (!endsWith(fdogDir, "/")) {
            fdogDir <- paste(fdogDir, "/", sep = "")
        }
        outPath <- paste(fdogDir, genomeName, sep = "")
    } else {
        outPath <- paste(
            root, "fdogout", "/", genomeName,
            sep = ""
        )
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
    
    if (reFdog == TRUE) {
        for (
            folder in list.dirs(
                outPath, recursive = TRUE, full.names = TRUE
            )
        ) {
            unlink(folder, recursive = TRUE)
        }
    }
    
    fastaSet <- lapply(
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
                            ".extended.fa", sep = ""
                        )
                    )
                ) {
                    exFasta <- readLines(
                        paste(
                            outPath, "/", coreGene, "/", coreGene, 
                            ".extended.fa", sep = ""
                        )
                    )
                    fastaList <- updateExFasta(
                        root, coreSet, exFasta, coreGene, genomeName, refSpec,
                        scoreMode
                    )
                    return(fastaList)
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
                        outPath, "/", coreGene, "/", coreGene, ".extended.fa", 
                        sep = ""
                    )
                )
            ){
                return(NULL)
            } else {
                exFasta <- readLines(
                    paste(
                        outPath, "/", coreGene, "/", coreGene, ".extended.fa",
                        sep = ""
                    )
                )
                fastaList <- updateExFasta(
                    root, coreSet, exFasta, coreGene, genomeName, refSpec,
                    scoreMode
                )
                return(fastaList)
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
    
    orthoMax <- 1
    for (fastaList in fastaSet) {
        if (length(fastaList) > orthoMax) {
            orthoMax <- length(fastaList)
        }
    }
    
    i <- 1:orthoMax
    splitedFasta <- lapply(
        i,
        function(i, fastaSet) {
            subFasta <- c()
            for (fastaList in fastaSet) {
                if (length(fastaList) >= i) {
                    subFasta <- c(subFasta, fastaList[[i]])
                }
            }
            return(subFasta)
        },
        fastaSet = fastaSet
    )
    
    if (!is.null(ppDir)) {
        outDir <- ppDir
    } else {
        outDir <- paste(
            root, "output", "/",  setName, "/", as.character(scoreMode), "/", 
            sep = ""
        )
    }
    
    temporary <- paste(outDir, "temporary", sep = "")
    if (!dir.exists(temporary)) {
        dir.create(temporary, recursive = TRUE)
    }
    
    o <- 1
    for (fasta in splitedFasta) {
        writeLines(
            fasta, 
            paste(temporary, "/", genomeName, "_", o, ".extended.fa", sep = "")
        )
        
        command <- paste(
            "fdogFAS",
            "-i", paste(
                temporary, "/", genomeName, "_", o, ".extended.fa", sep = ""
            ),
            "-w", weightPath,
            "-n", paste(genomeName, "_", o, sep = ""),
            "-o", temporary,
            "--cores", cpu
        )
        system(command)
        o <- o + 1
    }
    
    if (cleanup == TRUE) {
        unlink(outPath, recursive = TRUE)
    }
    
    outputList <- concanateFiles(temporary, genomeName)
    pp <- outputList[[1]]
    
    if (extend == TRUE) {
        domain0 <- outputList[[3]]
        domain1 <- outputList[[4]]
        exFasta <- outputList[[2]]
        
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
                paste(outDir, setName, "_forward.domains", sep = "")
            )
        ) {
            oriDomain0 <- read.table(
                paste(outDir, "_forward.domains", sep = ""),
                sep = "\t",
                comment.char = ""
            )
            oriDomain1 <- rbind(oriDomain1, domain1)
        } else {
            oriDomain1 <- domain1
        }
        
        if (
            file.exists(
                paste(outDir, setName, "_reverse.domains", sep = "")
            )
        ) {
            oriDomain1 <- read.table(
                paste(outDir, "_reverse.domains", sep = ""),
                sep = "\t",
                comment.char = ""
            )
            oriDomain0 <- rbind(oriDomain0, domain0)
        } else {
            oriDomain0 <- domain0
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
        write.table(
            oriDomain1,
            paste(outDir, setName, "_forward.domains", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
        write.table(
            oriDomain0,
            paste(outDir, setName, "_reverse.domains", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
        writeLines(
            oriFasta,
            paste(outDir, setName, ".extended.fa", sep = "")
        )
    }
    unlink(temporary, recursive = TRUE)
    return(pp)
}
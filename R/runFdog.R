#' fDOG in runFdog or runFdogBusco will be run with option --fasoff. The
#' extended fasta file of each core group will be merged and the merged file 
#' will be used as the input for fdogFAS. But the raw output from fDOG can not 
#' be used direct for the merging and must be processed. This function process
#' the extended fasta file of each core group depent on the using score mode.
#' For score mode 2 and 3 this function take the ortholog sequence from the 
#' extended fasta file and save it into a vector, the sequence of the references
#' species will be appended followed into the vector. If it exists more than one
#' ortholog sequence, the function will save the second ortholog sequence within
#' the references sequence into the other vector. For score mode 1 is analog but
#' instead of appending the references sequence into the vector, the function 
#' will append all training sequence of the core group into the vector. The 
#' function will returns a list, which contains the vector of the sequences. 
#' The number of the element of the list equal to the number of the orthologs,
#' that fDOG founded for this core group. This function will be used as a modul
#' in the function runFdog (not runFdogBusco).
#' 
#' @param root The path to the core directory, where the core set is stored 
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can 
#' contains more than one core set and the user must specify the interested 
#' core set. The core set will be stored in the folder core_orthologs in 
#' subfolder, specify them by the name of the subfolder
#' @param fasta The extended fasta file of the core group in form of a vector
#' @param coreGene The core group ID in the core set
#' @param genomeName The genome ID of the interested genome, Exp:HUMAN@9696@3
#' @param refSpec The genome ID of the references species
#' @param scoreMode the mode determines the method to scoring the founded 
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' 
#' @return A list, which contains the vector of the sequences. 
#' The number of the element of the list equal to the number of the orthologs,
#' that fDOG founded for this core group
#' @export

updateExFasta <- function(
    root, coreSet, fasta, coreGene, genomeName, refSpec, scoreMode
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

#' The output of fdogFAS is a phylogenetic profile (pp). But this pp still can
#' not be used to assess the completeneness of the interested genome. The pp
#' will contains the ortholog sequence within the training sequence. The FAS
#' score of each sequence is the FAS score between itself and the core group's
#' ortholog sequence. Therefore the score must be recalculated, for score mode
#' 2 and 3 the score between the ortholog sequence and the score of the training
#' sequence must be swaped. For score mode 1 the avarage score of the training 
#' sequences must be calculated and updated into the line of the ortholog
#' sequence. This function will be used as a modul in concanateFiles. The 
#' function will remove all the line of the training sequences and return the pp
#' of the interested genome.
#' 
#' @param pp The phylogenetic profile output from fdogFAS in form of a 
#' data.frame
#' @param genomeName The genome ID of the interested genome
#' 
#' @return The phylogenetic profile of the interested genome.
#' @export
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

#' After fdogFAS is finished. It will leave in a temporary folder the 
#' phylogenetic profiles (pp). The pps can not be used direct to assess the
#' completeness of the interested genome. This function processes the pp files
#' in a folder and merged them.
#'
#' @param directory path to the folder, that contains the pp files
#' @param genomeName the genomeID of the interested genome
#' @param scoreMode the mode determines the method to scoring the founded 
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @return a list which contains the pp of the interested genome in data.frame, 
#' the extended fasta file of the interested genome in form of a vector, which
#' contains the lines of the file and 2 domains file in form of data.frame
#' @export
concanateFiles <- function(directory, genomeName, scoreMode) {
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
            if (scoreMode == 1) {
                singlePP <- recalculateScore(singlePP, genomeName)
            }
            if (is.null(pp)) {
                pp <- singlePP
            } else {
                pp <- rbind(pp, singlePP)
            }
        }
    }
    pp <- extractPP(pp, genomeName)
    exFasta <- extractFasta(exFasta, genomeName)
    if (scoreMode != 1) {
        domain0 <- extractDomains(domain0, genomeName)
        domain1 <- extractDomains(domain1, genomeName)
    }
    return(list(pp, exFasta, domain0, domain1))
}

#' This function is used as a modul in the function computeReport() when using
#' socre mode 1, 2, 3. The function computeReport() take the path to the genome 
#' fasta file and the annotation file (if the user did not input it, it will be 
#' computed) of the interested genome as one of its inputs and arrange the 
#' files into the equivalent folder of the core directory. The genome fasta 
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
                    outPath, "/", refSpec, "/", coreGene, "/", 
                    "hamstrsearch.log", sep = ""
                )
            )) {
                if (
                    file.exists(
                        paste(
                            outPath, "/", refSpec, "/", coreGene, "/", coreGene, 
                            ".extended.fa", sep = ""
                        )
                    )
                ) {
                    if (scoreMode == 1) {
                        exFasta <- readLines(
                            paste(
                                outPath, "/", refSpec, "/", coreGene, "/", 
                                coreGene, ".extended.fa", sep = ""
                            )
                        )
                        fastaList <- updateExFasta(
                            root, coreSet, exFasta, coreGene, genomeName, 
                            refSpec, scoreMode
                        )
                        return(fastaList)
                    } else {
                        return(NULL)
                    }
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
                        ".extended.fa", sep = ""
                    )
                )
            ){
                return(NULL)
            } else {
                if (scoreMode == 1) {
                    exFasta <- readLines(
                        paste(
                            outPath, "/", refSpec, "/", coreGene, "/", coreGene, 
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
    
    if (scoreMode == 1) {
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
    }
    
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
    
    if (scoreMode == 1) {
        o <- 1
        for (fasta in splitedFasta) {
            writeLines(
                fasta, 
                paste(
                    temporary, "/", genomeName, "_", o, ".extended.fa", sep = ""
                )
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
    } else {
        o <- 1
        for (
            subfolder in list.dirs(
                outPath, full.names = TRUE, recursive = FALSE
            )
        ) {
            exFasFiles <- list.files(
                subfolder,
                pattern = "\\.extended.fa$",
                recursive = TRUE,
                full.names = TRUE
            )
            subcommand <- paste(exFasFiles, collapse = " ")
            command <- paste(
                "cat",
                subcommand,
                ">",
                paste(
                    temporary, "/", genomeName, "_", o, ".extended.fa", sep = ""
                )
            )
            system(command)
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
    }
    
    if (cleanup == TRUE) {
        unlink(outPath, recursive = TRUE)
    }
    
    outputList <- concanateFiles(temporary, genomeName, scoreMode)
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
            oriDomain1 <- read.table(
                paste(outDir, setName, "_forward.domains", sep = ""),
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
            oriDomain0 <- read.table(
                paste(outDir, setName, "_reverse.domains", sep = ""),
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
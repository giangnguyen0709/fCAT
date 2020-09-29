#' The function correct information in a original domain file. This function is
#' just applied when scoreMode=1
#'
#' @param domain The path to the domain file
#'
#' @return none
#' @export
extendDomain <- function(domain) {
    domainFile <- read.table(
        paste(domain, sep = ""), sep = "\t", comment.char = ""
    )
    compareSet <- unique(as.vector(domainFile$V1))
    domainSet <- lapply(compareSet,
        function(compare, domainFile) {
            sub <- subset(domainFile, V1 == compare)
            splited <- strsplit(compare, "#", fixed = TRUE)[[1]]
            splited2 <- strsplit(splited[1], "|", fixed = TRUE)[[1]]
            sub$V1 <- paste(splited2[1], "#", splited[2], sep = "")
            return(sub)
        },
        domainFile = domainFile
    )
    domainFile <- do.call("rbind", domainSet)
    write.table(domainFile, domain,
        row.names = FALSE,
        sep = "\t", quote = FALSE, col.names = FALSE
    )
}

#' The function calculate the FAS score in the scoreMode 1 and update the scores
#' into phylogenetic profile of a specific core gene
#'
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene
#' @param extend the logical argument to decide if the information of the
#' interested genome will be appended to the original files
#' @param refSpec the references species
#'
#' @return none
#' @export
getFasScore <- function(root, coreSet, coreGene, extend = FALSE, refSpec) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    setName <- coreSet
    queryGenome <- list.dirs(paste(root, "genome_dir", sep = ""),
        recursive = FALSE,
        full.names = FALSE
    )[[1]]

    seedFasta <- paste(root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
        coreGene, ".fa",
        sep = ""
    )
    queryFasta <- paste(root, "phyloprofile", "/", coreSet, "/", "1", "/",
        "hamstrout", "/", coreGene, "/", coreGene,
        ".extended.fa",
        sep = ""
    )
    lines <- readLines(queryFasta)
    newLines <- extractFasta(queryFasta, queryGenome)
    writeLines(newLines, queryFasta)

    seedFastaLines <- readLines(seedFasta)
    i <- 1:length(seedFastaLines)
    seedSet <- lapply(
        seedFastaLines[i[i %% 2 == 1]],
        function(seedID) {
            return(substr(seedID, 2, nchar(seedID)))
        }
    )
    seedSet <- unlist(seedSet)
    refID <- NULL
    for (seed in seedSet) {
        compare <- strsplit(seed, "|", fixed = TRUE)[[1]][2]
        if (compare == refSpec) {
            refID <- seed
            break
        }
    }
    annoDir <- paste(
        root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
        "fas_dir", "/", "annotation_dir",
        sep = ""
    )

    queryPP <- read.table(
        paste(
            root, "phyloprofile", "/", coreSet, "/", "1", "/", 
            "hamstrout", "/", coreGene, "/", coreGene, ".phyloprofile",
            sep = ""
        ), 
        header = TRUE, 
        sep = "\t"
    )
    queryPP <- extractPP(queryPP, queryGenome)
    FAS_F <- list()
    FAS_B <- list()

    for (orthoID in queryPP$orthoID) {
        FAS_F[as.character(orthoID)] <- 0
        FAS_B[as.character(orthoID)] <- 0
    }

    orthologNumber <- length(queryPP$orthoID)
    queryGenome <- strsplit(
        as.character(queryPP$orthoID[1]),
        "|",
        fixed = TRUE
    )[[1]][2]


    jobname <- paste(
        root, "phyloprofile", "/", coreSet, "/", "1", "/",
        "hamstrout", "/", coreGene,
        sep = ""
    )

    R.utils::createLink(
        paste(annoDir, "/", queryGenome, ".json", sep = ""),
        paste(root, "weight_dir", "/", queryGenome, ".json", sep = ""),
        overwrite = TRUE
    )
    # calculate FAS scores
    preScores <- lapply(seedSet,
        function(
            seedID, seedFasta, queryFasta, annoDir, root, coreSet, setName, 
            refID, jobname, queryGenome, coreGene
        ) {
            splited <- strsplit(seedID, "|", fixed = TRUE)[[1]]
            R.utils::createLink(
                paste(annoDir, "/", splited[2], ".json", sep = ""),
                paste(root, "weight_dir", "/", splited[2], ".json", sep = ""),
                overwrite = TRUE
            )
            command <- paste(
                "calcFAS",
                "-q", queryFasta,
                "-s", seedFasta,
                "-o", jobname,
                "-a", annoDir,
                "--seed_id", paste('"', seedID, '"', sep = ""),
                "-t", 10,
                "--ref_2", paste(root, "genome_dir",
                    "/", splited[2], "/",
                    splited[2], ".fa",
                    sep = ""
                ),
                "-r", paste(root, "genome_dir",
                    "/", queryGenome, "/",
                    queryGenome, ".fa",
                    sep = ""
                ),
                "--bidirectional",
                "-n", coreGene,
                "--raw"
            )
            if (seedID != refID) {
                command <- paste(command, "--tsv")
            }
            lines <- system(command, intern = TRUE)
            file.remove(paste(annoDir, "/", splited[2],
                ".json",
                sep = ""
            ))
            return(lines)
        },
        queryFasta = queryFasta,
        seedFasta = seedFasta,
        root,
        coreSet = coreSet,
        setName = setName,
        refID = refID,
        jobname = jobname,
        queryGenome = queryGenome,
        annoDir = annoDir,
        coreGene = coreGene
    )
    file.remove(paste(annoDir, "/", queryGenome, ".json", sep = ""))
    if (extend == TRUE) {
        try(extendDomain(
            paste(jobname, "/", coreGene, "_forward.domains", sep = "")),
            silent = TRUE
        )
        try(extendDomain(
            paste(jobname, "/", coreGene, "_reverse.domains", sep = "")),
            silent = TRUE
        ) 
    }

    bcount <- 0
    fcount <- 0
    for (lines in preScores) {
        count <- 1
        for (line in lines) {
            if (startsWith(line, "#")) {
                splited <- strsplit(line, "\t", fixed = TRUE)[[1]]
                if (count > orthologNumber) {
                    q <- strsplit(splited[2], "|", fixed = TRUE)[[1]][2]
                    s <- strsplit(splited[3], "|", fixed = TRUE)[[1]][2]
                    score <- as.numeric(splited[length(splited)])
                    FAS_B[[q]] <- FAS_B[[q]] + score
                    bcount <- bcount + 1
                } else {
                    q <- strsplit(splited[3], "|", fixed = TRUE)[[1]][2]
                    s <- strsplit(splited[2], "|", fixed = TRUE)[[1]][2]
                    score <- as.numeric(splited[length(splited)])
                    FAS_F[[q]] <- FAS_F[[q]] + score
                    fcount <- fcount + 1
                }
            }
        }
    }
    FAS_F <- unlist(FAS_F)
    FAS_B <- unlist(FAS_B)
    FAS_F <- FAS_F / fcount
    FAS_B <- FAS_B / bcount
    pp <- cbind(queryPP, FAS_F, FAS_B)
    write.table(pp, 
                paste(
                    root, "phyloprofile", "/", coreSet, "/", "1", "/", 
                    "hamstrout", "/", coreGene, "/", coreGene, ".phyloprofile",
                    sep = ""),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
    )
}

#' The function calculate the FAS score in the scoreMode 1 and update the scores
#' into phylogenetic profile of a specific core gene
#'
#' @param root Path to the root folder
#' @param coreSet The core set name
#' @param coreGene the ID of the core gene
#' @param extend the logical argument to decide if the information of the
#' interested genome will be appended to the original files
#' @param refSpec the references species
#'
#' @return none
#' @export
updateFasScore <- function(root, coreSet, coreGene, extend, refSpec) {
    if (!endsWith(coreSet, "/")) {
        coreSet <- paste(coreSet, "/", sep = "")
    }
    getFasScore(root, coreSet, coreGene, extend, refSpec)
    file.remove(
        paste(
            root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
            "fas_dir", "/", "annotation_dir", "/", coreGene, "extended.json",
            sep = ""
        )
    )
    unlink(
        paste(
            root, "core_orthologs", "/", coreSet, "/", coreGene, "/", "fas_dir", 
            "/", "annotation_dir", "/", "tmp", 
            sep = ""
        ),
        recursive = TRUE
    )
}

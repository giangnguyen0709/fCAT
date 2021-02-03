#' The function creates FAS annotation for the fasta file of a core gene 
#' in the core set. This annotation file will be saved in each core gene folder
#' and be used then to calculate the cut off value for each score mode.
#'
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param coreGene The ID of the core gene in the set
#' @return none
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' createAnnotation(coreFolder, "test", "530670")
#' @export
createAnnotation <- function(root, coreSet, coreGene) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }

    annoFolder <- paste(
        root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
        "fas_dir", "/", "annotation_dir",
        sep = ""
    )
    if (!dir.exists(annoFolder)) {
        dir.create(annoFolder, recursive = TRUE)
    }

    command <- paste(
        "annoFAS",
        "-i", paste(root, "core_orthologs", "/", coreSet, "/",
            coreGene, "/", coreGene, ".fa",
            sep = ""
        ),
        "-o", annoFolder,
        "-n", paste(coreGene, sep = "")
    )
    system(command)
    tmp <- paste(annoFolder, "/", "tmp", sep = "")
    if (dir.exists(tmp)) {
        unlink(tmp, recursive = TRUE)
    }
}

#' The function run createAnnotation function in a loop to compute the FAS
#' annotation file for all core genes in the core set
#'
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#'
#' @return none
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' createAllAnnotation(coreFolder, "test")
#' @export
createAllAnnotation <- function(root, coreSet) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    lapply(
        list.dirs(paste(root, "core_orthologs", "/", coreSet, sep = ""),
        recursive = FALSE, full.names = FALSE
    ),
    function(coreGene, root, coreSet) {
        print(paste("Annotating for", coreGene, sep = " "))
        createAnnotation(root, coreSet, coreGene)
    },
    root = root,
    coreSet = coreSet
    )
}

#' The function calculate all cutoff values for all mode of a specific core 
#' gene in the set. For score mode 1 it will calculate the avarage of all vs 
#' all FAS scores between the training sequences in the core gene. For score 
#' mode 2 it will calculate the avarage of the FAS score between each sequence 
#' against all training sequences in the core gene. The scores will be writen 
#' in a table with a column is the ID of the sequences and a column is the 
#' corresponding value. For score mode 3, the function will calculate the 
#' avarage of 1 vs all FAS scores for each training sequence in the core gene. 
#' The avarages build a distribution, the function will calculate the 
#' confidence interval of this distribution and write the upper value and the 
#' lower value of the interval in a file in the core gene folder.
#'
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param coreGene the ID of the core gene in the set
#'
#' @return none
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' calculateCutoff(coreFolder, "test", "530670")
#' @export
calculateCutoff <- function(root, coreSet, coreGene) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    fastaFile <- paste(
        root, "core_orthologs", "/", coreSet, "/",
        coreGene, "/", coreGene, ".fa",
        sep = ""
    )
    annoDir <- paste(
        root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
        "fas_dir", "/", "annotation_dir",
        sep = ""
    )

    fasta <- readLines(fastaFile)
    i <- 1:length(fasta)
    i <- i[i %% 2 == 1]
    querySet <- lapply(
        fasta[i],
        function(header) {
            return(substr(header, 2, nchar(header)))
        }
    )

    genomeSet <- lapply(querySet, function(query) {
        return(strsplit(query, "|", fixed = TRUE)[[1]][2])
    })

    scoreDist <- list()
    genomeScores <- list()
    for (i in 1:length(genomeSet)) {
        genomeScores[[genomeSet[[i]]]] <- 0
        if (i != length(genomeSet)) {
            scoreDist[[genomeSet[[i]]]] <- list()
        }
    }

    for (i in 1:(length(genomeSet) - 1)) {
        for (j in (i + 1):length(genomeSet)) {
            scoreDist[[genomeSet[[i]]]][[genomeSet[[j]]]] <- 0
        }
    }

    scoreSet <- lapply(
        querySet,
        function(queryID, fastaFile, annoDir, root, coreSet) {
            refSpec <- strsplit(queryID, "|", fixed = TRUE)[[1]][2]
            refProteome <- paste(
                root, "weight_dir", "/",
                refSpec, ".json",
                sep = ""
            )
            R.utils::createLink(
                paste(annoDir, "/", refSpec, ".json", sep = ""),
                refProteome,
                overwrite = TRUE
            )
            refProteome <- paste(
                root, "genome_dir", "/", refSpec, "/", refSpec, ".fa", sep = "")
            command <- paste(
                "calcFAS",
                "-q", fastaFile,
                "-s", fastaFile,
                "--query_id", paste('"', queryID,
                    '"',
                    sep = ""
                ),
                "-a", annoDir,
                "-o", annoDir,
                "--tsv",
                "-r", refProteome,
                "-t", "10",
                "--raw"
            )
            
            lines <- system(command, intern = TRUE)

            scores <- list()
            for (line in lines) {
                if (startsWith(line, "#")) {
                    splited <- strsplit(line, "\t", fixed = TRUE)[[1]]
                    q <- strsplit(splited[3], "|", fixed = TRUE)[[1]][2]
                    s <- strsplit(splited[2], "|", fixed = TRUE)[[1]][2]
                    score <- as.numeric(splited[length(splited)])
                    pack <- list(q, s, score)
                    scores[[length(scores) + 1]] <- pack
                }
            }
            file.remove(paste(annoDir, "/", refSpec,
                ".json",
                sep = ""
            ))
            return(scores)
        },
        fastaFile,
        annoDir,
        root = root,
        coreSet
    )
    
    for (level1 in scoreSet) {
        for (score in level1) {
            if (score[[1]] != score[[2]]) {
                a <- genomeScores[[score[[1]]]]
                b <- genomeScores[[score[[2]]]]
                genomeScores[[score[[1]]]] <- a + score[[3]]
                genomeScores[[score[[2]]]] <- b + score[[3]]
                if (!is.null(scoreDist[[score[[1]]]][[score[[2]]]])) {
                    v <- scoreDist[[score[[1]]]][[score[[2]]]]
                    scoreDist[[score[[1]]]][[score[[2]]]] <- v + score[[3]]
                } else {
                    v <- scoreDist[[score[[2]]]][[score[[1]]]]
                    scoreDist[[score[[2]]]][[score[[1]]]] <- v + score[[3]]
                }
            }
        }
    }
    scoreDist <- unlist(scoreDist) / 2
    genomeScores <- unlist(genomeScores) / ((length(genomeSet) - 1) * 2)
    
    avaMean <- mean(genomeScores)
    
    features <- EnvStats::eexp(scoreDist, ci = TRUE)
    lcl <- 1 / (features$interval$limits[[2]])
    ucl <- 1 / (features$interval$limits[[1]])

    scoreFolder <- paste(
        root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
        "fas_dir", "/", "score_dir",
        sep = ""
    )
    if (!dir.exists(scoreFolder)) {
        dir.create(scoreFolder, recursive = TRUE)
    }
    # Table 1:
    label <- c("mean", "LCL", "UCL")
    value <- c(avaMean, lcl, ucl)
    cutoffTable <- data.frame(label, value)

    cutoffFile <- paste(scoreFolder, "/", "1", ".cutoff", sep = "")
    write.table(
        cutoffTable, cutoffFile,
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    # Table 2:
    cutoffTable <- data.frame(
        taxa = unlist(genomeSet),
        cutoff = genomeScores
    )
    cutoffFile <- paste(scoreFolder, "/", "2", ".cutoff", sep = "")
    write.table(
        cutoffTable, cutoffFile,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}

#' The function run calculateCutoff in a loop to calculate cutoff 
#' values for all core genes in the core set
#'
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#'
#' @return none
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' calculateAllCutoff(coreFolder, "test")
#' @export
calculateAllCutoff <- function(root, coreSet) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }

    coreOrtho <- paste(root, "core_orthologs", "/", coreSet, "/", sep = "")
    lapply(list.dirs(coreOrtho, full.names = FALSE, recursive = FALSE),
        function(coreGene, root, coreSet, mode) {
            print(paste("Starting calculate cutoff for", coreGene, sep = " "))
            calculateCutoff(root, coreSet, coreGene)
        },
        root = root,
        coreSet = coreSet
    )
}

#' The function run createAllAnnotation and calculateAllCufoff to process
#' the core set
#'
#' @param coreDir The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#'
#' @return none
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' processCoreSet(coreFolder, "test")
#' @export
processCoreSet <- function(coreDir, coreSet) {
    startTime <- Sys.time()
    if (!endsWith(coreDir, "/")) {
        coreDir <- paste(coreDir, "/", sep = "")
    }

    createAllAnnotation(coreDir, coreSet)
    calculateAllCutoff(coreDir, coreSet)
    endTime <- Sys.time()
    print("Running time:")
    print(endTime - startTime)
}

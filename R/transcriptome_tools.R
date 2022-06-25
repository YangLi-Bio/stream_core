# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
data.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/"


# Source the R scripts
source(paste0(code.dir, "transcriptome_tools.R"))


# LTMG modeling
#' @import IRISFGM
call_LTMG <- function(obj) {

  # Libraries
  library(IRISFGM)


  ltmg.df <- as.data.frame(obj@assays$RNA@data) # get the expression matrix saved in a data frame
  ltmg.obj <- ProcessData(CreateIRISFGMObject(ltmg.df),
                          normalization = "cpm", IsImputation = F)
  ltmg.matrix <- GetLTMGmatrix(RunLTMG(ltmg.obj, Gene_use = "all"))
  ltmg.matrix <- ltmg.matrix - min(ltmg.matrix) # avoid the overflow of 1s
  message ("Finished LTMG modeling for ", nrow(ltmg.matrix), " variable genes and ",
           ncol(ltmg.matrix), " cells.\n")


  return(ltmg.matrix)
}


# Write LTMG matrix
write_LTMG <- function(LTMG.matrix, p.LTMG) {

  writeLines("o", con = p.LTMG, sep = "\t") # print the beginning symbol, i.e., "o"
  write.table(LTMG.matrix, p.LTMG, append = T, quote = F, sep = "\t") # write the matrix
  message ("Finished writing the discretized gene expression to file ", p.LTMG, ".\n")
}


# # Run QUBIC for biclustering
# run_QUBIC <- function(path = "/fs/ess/PCON0022/liyang/tools/biclustering/QUBIC2-master",
#                       file, c = 1.0, dual = F, k = 50, o = 100, f = 0.25,
#                       d = T, square = T) {
#
#   command.ln <- paste0 (path, "/qubic -i ",
#                         file, " -c ", c, " -k ", k, " -o ", o, " -f ", f)
#   message ("--------> QUBIC2 command: ", command.ln, "\n")
#
#   if (d) { # performed on discretized matrix
#     command.ln <- paste0 (command.ln, " -d")
#   }
#
#   if (dual) { # perform dual expansion
#     command.ln <- paste0 (command.ln, " -C")
#   }
#
#   if (square) { # perform dual expansion
#     command.ln <- paste0(command.ln, " -N")
#   }
#
#   system (command.ln, wait = T)
#
#   return (1)
# }


# Select the next bicluster
#' @import pbmcapply
select_block <- function(block.original, retained.blocks, ucells = c(),
                         ugenes = NULL, sim.mode = "both") {

  # Libraries
  library(pbmcapply)

  add.size.ll <- pbmclapply(seq_along(block.original), function(i) {
    if (i %in% retained.blocks) {
      return(-1)
    }

    ifelse (sim.mode == "both",
            add.size <- length(setdiff(block.original[[i]]$genes, ugenes)) *
              length(setdiff(block.original[[i]]$cells, ucells)),
            add.size <- length(setdiff(block.original[[i]]$cells, ucells)))
    # compute the add size

    return(add.size)
  }, mc.cores = min(detectCores(), length(block.original))) %>% unlist

  max.add.value <- max(add.size.ll)
  max.add.id <- which.max(add.size.ll)

  if (max.add.value == -1) {
    return(-1)
  } else {
    return(max.add.id)
  }
}


# Retain QUBIC blocks based on submodular optimization
retain_blocks <- function(block.original, sim.mode = "both",
              cover.blocks = 10) {
  i <- 0
  retained.blocks <- c() # selected blocks
  ucells <- c() # an empty set
  ugenes <- NULL # an empty set

  if (sim.mode == "both") {
    ugenes <- c()
  }

  while(i < cover.blocks) {
    block.id <- select_block(block.original = block.original,
                             retained.blocks = retained.blocks,
                             ucells = ucells,
                             ugenes = ugenes,
                             sim.mode = sim.mode)
    # select the next block

    if (block.id == -1) {
      break
    }

    retained.blocks <- c(retained.blocks, block.id)
    i <- i + 1
    ucells <- union(ucells, block.original[[block.id]]$cells)

    if (sim.mode == "both") {
      ugenes <- union(ugenes, block.original[[block.id]]$genes)
    }
  }

  retained.blocks
}

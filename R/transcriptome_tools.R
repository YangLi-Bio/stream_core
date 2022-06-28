# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
data.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/"


# Source the R scripts


# LTMG modeling
#' @import IRISFGM
call_LTMG <- function(obj) {

  # Libraries
  library(IRISFGM)


  ltmg.df <- as.data.frame(obj@assays$RNA@data) # get the expression matrix saved in a data frame
  ltmg.obj <- ProcessData(CreateIRISFGMObject(ltmg.df),
                          normalization = "cpm", IsImputation = F)
  # ltmg.matrix <- GetLTMGmatrix(RunLTMG(ltmg.obj, Gene_use = "all"))
  # ltmg.matrix <- ltmg.matrix - min(ltmg.matrix) # avoid the overflow of 1s
  # message ("Finished LTMG modeling for ", nrow(ltmg.matrix), " variable genes and ",
  #          ncol(ltmg.matrix), " cells.\n")
  #
  #
  # return(ltmg.matrix)


  RunLTMG(ltmg.obj, Gene_use = "all")
}


# Write LTMG matrix
write_LTMG <- function(LTMG.matrix, p.LTMG) {

  writeLines("o", con = p.LTMG, sep = "\t") # print the beginning symbol, i.e., "o"
  write.table(LTMG.matrix, p.LTMG, append = T, quote = F, sep = "\t") # write the matrix
  message ("Finished writing the discretized gene expression to file ", p.LTMG, ".\n")
}


# Run QUBIC for biclustering
run_QUBIC <- function(path, file, c = 1.0, dual = F, k = 50, o = 100, f = 0.25,
                      d = T, square = T) {

  command.ln <- paste0 (path, " -i ",
                        file, " -c ", c, " -k ", k, " -o ", o, " -f ", f)
  message ("QUBIC2 command: ", command.ln, "\n")
  if (d) { # performed on discretized matrix
    command.ln <- paste0 (command.ln, " -d")
  }
  if (dual) { # perform dual expansion
    command.ln <- paste0 (command.ln, " -C")
  }
  if (square) { # perform dual expansion
    command.ln <- paste0(command.ln, " -N")
  }
  system (command.ln, wait = T)


  return (1)
}


# Select the next bicluster
#' @import pbmcapply
select_block <- function(block.list, retained.blocks, ucells = c(),
                         ugenes = NULL, sim.mode = "both") {

  # Libraries
  library(pbmcapply)

  add.size.ll <- pbmclapply(seq_along(block.list), function(i) {
    if (i %in% retained.blocks) {
      return(-1)
    }

    ifelse (sim.mode == "both",
            add.size <- length(setdiff(block.list[[i]]$genes, ugenes)) *
              length(setdiff(block.list[[i]]$cells, ucells)),
            add.size <- length(setdiff(block.list[[i]]$cells, ucells)))
    # compute the add size

    return(add.size)
  }, mc.cores = min(detectCores(), length(block.list))) %>% unlist

  max.add.value <- max(add.size.ll)
  max.add.id <- which.max(add.size.ll)

  if (max.add.value == -1) {
    return(-1)
  } else {
    return(max.add.id)
  }
}


# Retain QUBIC blocks based on submodular optimization
retain_blocks <- function(block.list, sim.mode = "both",
              cover.blocks = 10) {
  i <- 0
  retained.blocks <- c() # selected blocks
  ucells <- c() # an empty set
  ugenes <- NULL # an empty set

  if (sim.mode == "both") {
    ugenes <- c()
  }

  while(i < cover.blocks) {
    block.id <- select_block(block.list = block.list,
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
    ucells <- union(ucells, block.list[[block.id]]$cells)

    if (sim.mode == "both") {
      ugenes <- union(ugenes, block.list[[block.id]]$genes)
    }
  }

  retained.blocks
}


# Get the list of QUBIC biclusters
#' @import dplyr
get_block_list <- function(block.genes, block.cells, ngene = 3,
                           top.block = NULL, top.cells = NULL) {

  # Libraries
  library(dplyr)


  if (is.null(top.block)) {
    top.block <- length(unique(block.genes$Label))
  }
  gene.ll <- split(x = block.genes, f = block.genes$Label) %>%
    lapply(., `[[`, ('Gene'))
  cell.ll <- split(x = block.cells, f = block.cells$Label) %>%
    lapply(., `[[`, ('Gene'))
  flag.ll <- sapply(gene.ll, length) %>% `>=` (ngene) # check the eligibility of each block
  list1 <- lapply(seq_along(gene.ll), function(i) {
    return(list(genes = gene.ll[[i]], cells = cell.ll[[i]]))
  })
  list3 <- list1[1 : min(top.block, length(list1))]


    if (is.null(top.cells)) {
    return(list3)
  }


  ncells.ll <- sapply(cell.ll, length)
  cell.order <- order(ncells.ll, decreasing = T)
  gene.ll <- gene.ll[cell.order]
  cell.ll <- cell.ll[cell.order]
  list2 <- lapply(seq_along(gene.ll), function(i) {
    return(list(genes = gene.ll[[i]], cells = cell.ll[[i]]))
  })
  list4 <- list2[1 : min(top.cells, length(list2))]


  return(c(list3, list4))
}


# load QUBIC blocks
get_block <-function(block.file, keyword = 'Genes', n = 200) {

  tmp.block <- readLines(block.file) # file.path is used to
  tmp.bc <- grep(keyword, tmp.block, value = T) # value: if FALSE, a vector containing
  tmp.cel.module <- sapply(strsplit(tmp.bc,':'),'[', 2) # '[' and '2' in sapply means outputting
  GENES <- as.character()   # store the genes
  label_G <- as.numeric()   # store the occurrence of genes
  for (j in 1:length(tmp.cel.module)) {
    if (j > n) {
      break
    }
    BCgene <- unlist(strsplit(tmp.cel.module[j], split = " ")) # transform a list into vector. Why?
    BCgene <- BCgene[BCgene != ""]  # exclude the blank string (maybe the last one?)
    GENES <- c(GENES, BCgene)
    label_G <- c(label_G, rep(j, length(BCgene)))
  }
  df_C <- data.frame(gene_name = GENES, label = as.factor(label_G))
  gene.list <- as.character(df_C$gene_name) # remove the double quotes?
  df_C$gene_name <- gene.list
  colnames(df_C) <- c("Gene", "Label")
  message (length(unique(df_C$Label)), ' QUBIC 2.0 blocks were identified.\n')


  return(df_C)
}


# Load QUBIC biclusters
load_blocks <- function(block.file = "./LTMG_matrix.txt",
            cover.blocks = 10, n.blocks = 200,
            sim.mode = "cell", rank.blocks = T) {

  block.genes <- get_block(block.file, n = n.blocks) # get block genes from QUBIC2 output
  block.cells <- get_block(block.file, keyword = 'Conds',
                           n = n.blocks) # get block cells from QUBIC2 output
  # block.genes <- block.original@BiCluster@CoReg_gene
  # block.cells <- block.original@BiCluster@CoCond_cell
  if (length(levels(block.genes$Label)) !=
      length(levels(block.cells$Label))) {
    stop ("The number of gene sets and cell sets in QUBIC 2.0 blocks are different: \n")
  }
  block.list <- get_block_list(block.genes = block.genes, block.cells = block.cells)
  if (length(block.list) < 1) {
    stop ("There is no QUBIC bicluster identified!\n")
  }

  if (rank.blocks) {
    return(block.list[retain_blocks(block.list = block.list, sim.mode = sim.mode,
                                  cover.blocks = cover.blocks)])
  } else {
    return(block.list[1 : min(cover.blocks, length(block.list))])
  }
}


# Write the LTMG matrix into txt file
writeLTMG <- function(LTMG.matrix = NULL, LTMG.file = "./") {
  writeLines("o", con = LTMG.file, sep = "\t") # print the beginning symbol, i.e., "o"
  write.table(LTMG.matrix, LTMG.file, append = T, quote = F, sep = "\t") # write the matrix
  message ("Wrote the discretized gene expression to file ", LTMG.file, ".\n")
}

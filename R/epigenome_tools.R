# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
data.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/"


# Link CREs to CREs via covariance
#' @import Signac GenomicRanges Matrix
link_cor <- function(x.atac, distance = 500000, cicero.covar = 0) {

  # Libraries
  library(Signac)
  library(GenomicRanges)
  library(Matrix)


  sep <- c("-", "-")
  GR.ll <- StringToGRanges(rownames(x.atac)) # transform strings to GRanges
  GR.expand <- resize(x = GR.ll, width = distance, fix = 'center') # find the TSS location
  overlaps <- findOverlaps(
    query = GR.ll,
    subject = GR.expand,
    type = 'any',
    select = 'all'
  ) # find the peaks overlaped with the extended genomic ranges of peaks
  message ("Finished finding overlapping peaks within ", distance, " bps.\n")
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = GR.ll), length(x = GR.expand))
  ) # build a sparse matrix to record the overlaps between peaks and extended genomic ranges of genes
  rownames(x = hit_matrix) <- GRangesToString(grange = GR.ll, sep = sep) # use peak names as the row names
  colnames(x = hit_matrix) <- GRangesToString(grange = GR.ll, sep = sep) # use peak names as the column names
  hit_matrix[!upper.tri(hit_matrix, diag = F)] <- 0 # mask the diagonal and upper triangular elements
  row.col <- which(hit_matrix == 1, arr.ind = T) # double vector to hold the row and column IDs for elements
  # equalling one
  if (nrow(row.col) < 1) {
    message ('No CRE-CRE linkage is identified.\n')
    return(NULL)
  }


  # Assign covariance as the matrix elements
  for (i in 1 : nrow(row.col)) {
    hit_matrix[row.col[i, 1], row.col[i, 2]] <- cov(x.atac[row.col[i, 1], ],
                                                    x.atac[row.col[i, 2], ])
  }
  hit_matrix[which(hit_matrix < cicero.covar)] <- 0 # discard the negative elements
  summ <- summary(hit_matrix) # convert the matrix into a sparse matrix
  cor.links <- data.frame(Origin = rownames(hit_matrix)[summ$i],
                          Destination = colnames(hit_matrix)[summ$j],
                          Weight      = summ$x) # transform the sparse matrix into a data frame
  colnames(cor.links) <- c('node1', 'node2', 'weight') # name the columns
  max.weight <- max(cor.links$weight)
  min.weight <- min(cor.links$weight)
  diff <- max.weight - min.weight
  quiet(ifelse(diff > 0, cor.links$weight <- (max.weight - cor.links$weight) /
                 diff, cor.links$weight <- 0)) # normalize the weights
  message ('Finished generating ', nrow(cor.links), ' CRE-CRE linkages.\n')


  cor.links
}


# Link CREs to CREs via Cicero
#' @import cicero dplyr BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg38 BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm9 data.table
link_cicero <- function(x, distance = 500000, cicero.covar = 0,
                        org.gs = BSgenome.Hsapiens.UCSC.hg38) {

  # Libraries
  library(cicero)
  library(dplyr)
  library(org.gs, character.only = T)
  library(data.table)


  summ <- summary(x) # convert the matrix into a sparse matrix
  cicero.data <- data.frame(Origin = rownames(x)[summ$i],
                            Destination = colnames(x)[summ$j],
                            Weight      = summ$x) # transform the sparse matrix into a data frame
  input.cds <- make_atac_cds(cicero.data, binarize = F) %>% detect_genes
  input.cds <- input.cds[Matrix::rowSums(exprs(input.cds)) != 0, ] %>% estimate_size_factors %>%
    preprocess_cds(method = "LSI", verbose = F) %>% reduce_dimension(reduction_method = 'UMAP',
                                                                     preprocess_method = "LSI")
  umap.coords <- reducedDims(input.cds)$UMAP # obtain the UMAP coordinates
  cicero.cds <- make_cicero_cds(input.cds, reduced_coordinates = umap.coords)
  genome.info <- data.frame(org.gs@seqinfo@seqnames,
                            org.gs@seqinfo@seqlengths) # genome sequence lengths
  colnames(genome.info) <- c("seqnames", "seqlengths") # rename the columns
  cicero.links <- run_cicero(cicero.cds, genome.info,
                             window = distance) # build peak-peak linkages using cicero
  colnames(cicero.links) <- c('node1', 'node2', 'weight')
  cicero.links$node2 <- as.character(cicero.links$node2) # convert factors into characters
  cicero.links <- rbindlist(apply(cicero.links, 1, function(r) {
    ifelse(r[1] >= r[2], return(list(node1 = r[1], node2 = r[2], weight = r[3])),
           return(list(node1 = r[2], node2 = r[1], weight = r[3])))
  }), fill = T) %>% dplyr::distinct()


  cicero.links$weight <- as.numeric(cicero.links$weight) # convert characters into numeric values
  cicero.links$weight[is.na(cicero.links$weight)] <- 0 # substitute the NA values
  max.weight <- max(cicero.links$weight)
  min.weight <- min(cicero.links$weight)
  diff <- max.weight - min.weight
  ifelse(diff > 0, cicero.links$weight <- (max.weight - cicero.links$weight) /
           diff, cicero.links$weight <- 0) # normalize the weights
  cat('Finished generating', nrow(cicero.links), 'CRE-CRE linkages.\n')

  cicero.links
}

# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"


# Source the R scripts


# Calculate the distance to the nearest TSS
#' @import GenomicRanges Matrix Signac
DistanceToTSS <- function(peaks, genes, distance = 200000,
                          sep = c("-", "-")) {

  # Libraries
  library(GenomicRanges)
  library(Matrix)


  tss <- resize(x = genes, width = 1, fix = 'start') # find the TSS location
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  ) # extand the genomic range from the TSS till downstream/upstream 200000 bp
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  ) # find the peaks overlaped with the extended genomic ranges of genes
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  ) # build a sparse matrix to record the overlaps between peaks and extended genomic ranges of genes
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep) # use peak names as the row names
  colnames(x = hit_matrix) <- genes.extended$gene_name # use gene names as the column names
  hit_matrix
}


# Retain the longest transcripts for all exons
#' @import data.table GenomicRanges
CollapseToLongestTranscript <- function(ranges) {

  # Libraries
  library(data.table)
  library(GenomicRanges)


  range.df <- as.data.table(x = ranges) # transform a GRanges object into a data frame
  range.df$strand <- as.character(x = range.df$strand) # transform GRanges into character strings
  range.df$strand <- ifelse( # check whether the strand information is available
    test = range.df$strand == "*", # ambiguous
    yes = "+", # treat all sequences by only using their positive strands
    no = range.df$strand # use the provided strands
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]]),
    "gene_name"
  ] # merge exons into the longest transcripts as a representative of the gene
  colnames(x = collapsed) <- c(
    "gene_name", "seqnames", "start", "end", "strand", "gene_biotype"
  )
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE # the information not used to build the data frame
    # will be retained in meta data
  )
  gene.ranges
}


#' @import Signac Seurat pbmcapply dplyr
#' @export
filter_nearby_genes <- function(obj, distance = 500000, peak.assay = "ATAC") {

  # Libraries
  library(Signac)
  library(Seurat)
  library(pbmcapply)
  library(dplyr)


  # calculate the nearby genes
  gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = obj[[peak.assay]]))
  peaks <- StringToGRanges(rownames(obj[[peak.assay]])) # get peak coordinates
  distance.df <- summary(DistanceToTSS(peaks = peaks, genes = gene.coords,
                                       distance = distance)) # distance matrix
  peak.names <- rownames(obj[[peak.assay]]) # get peak names
  gene.names <- gene.coords$gene_name # get gene names


  rbindlist(pbmclapply(1:nrow(distance.df), function(i) {
    return(list(peak = peak.names[distance.df[i, 1]], gene = gene.names[distance.df[i, 2]]))
  }, mc.cores = min(detectCores(), nrow(distance.df))), fill = T) %>%
    dplyr::filter(gene %in% rownames(obj[["RNA"]]))
}


# Link CREs to genes
#' @import Signac pbmcapply dplyr Seurat
link_signac <- function(x, distance = 500000,
                        signac.score = 0,
                        signac.pval = 1,
                        min.cells = 10,
                        peak.assay = 'ATAC') {

  # Libraries
  library(Seurat)
  library(Signac)
  library(pbmcapply)
  library(dplyr)


  DefaultAssay(x) <- peak.assay # set 'ATAC' as he default assay
  xx <- x # back up
  x <- tryCatch(LinkPeaks(object = x, distance = distance,
                           min.cells = min.cells,
                           peak.assay = peak.assay, expression.assay = 'RNA',
                           pvalue_cutoff = signac.pval,
                           score_cutoff = signac.score, verbose = T),
                error = function(e) {
                  0 }) # build linkages
  if (!is.numeric(x)) {
    signac.links <- data.frame(node1 = x[[peak.assay]]@links$peak,
                               node2 = x[[peak.assay]]@links$gene,
                               weight = x[[peak.assay]]@links$score)
  } else {
    x <- xx
    signac.links <- filter_nearby_genes(obj = x) # link peaks to genes using heuristics
    signac.links <- cbind(signac.links, pbmclapply(1 : nrow(signac.links),
                                                   function(i) {
      vx <- as.vector(x[["RNA"]][signac.links$gene[i]])
      vy <- as.vector(x[[peak.assay]][signac.links$peak[i]])

      if (sum(vx > 0) <= 0 | sum(vy > 0) <= 0) {
        return(0)
      } else if (sd(vx) == 0 | sd(vy) == 0) {
        return(1)
      } else {
        return(cor(x = vx, y = vy, method = "pearson"))
      }
    }, mc.cores = min(detectCores(), nrow(signac.links))) %>% unlist)
    colnames(signac.links) <- c("node1", "node2", "weight")
    signac.links <- signac.links[signac.links$weight > signac.score,] # filter linkages
  } # no peaks are linked to the genes
  if (nrow(signac.links) < 1) {
    return(NULL)
  }


  # Calculate weights
  max.weight <- max(signac.links$weight)
  min.weight <- min(signac.links$weight)
  diff <- max.weight - min.weight
  ifelse(diff > 0, signac.links$weight <- (max.weight - signac.links$weight) /
           diff, signac.links$weight <- 0) # normalize the weights
  message ('Finished generating ', nrow(signac.links), ' CRE-gene linkages.\n')

  signac.links
}


# Get the most coherent peak-gene pairs
#' @import data.table
get_coherent_peak_gene_pairs <- function(peak_distance_matrix,
                                         HBC.rna, HBC.atac) {

  # Libraries
  library(data.table)


  row.col <- which(peak_distance_matrix > 0, arr.ind = T) # get the nonzero elements
  peak.gene <- rbindlist(apply(row.col, 1, function(rr) {
    return(list(rownames(peak_distance_matrix)[rr[1]], colnames(peak_distance_matrix)[rr[2]]))
  })) # convert row and column ids into peaks and rows
  colnames(peak.gene) <- c("peak", "gene")
  peak.gene.weight <- cbind(peak.gene, weight = apply(peak.gene, 1, function(rr) {
    return(length(which(HBC.atac[rr[1], ] > 0 & HBC.rna[rr[2], ] > 0)))
  }))
  peak.gene.weight <- peak.gene.weight[with(peak.gene.weight, order(gene, -weight))] # sort the data frame


  return(tapply(peak.gene.weight$peak, peak.gene.weight$gene, function(i) {
    head(i, n = 1) # select the pairs of peaks and genes
  }))
}

# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"


# Source the r scripts


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


#' @import Signac Seurat pbmcapply
#' @export
filter_nearby_genes <- function(obj, distance = 500000, peak.assay = "ATAC") {
  
  # Libraries
  library(Signac)
  library(Seurat)
  library(pbmcapply)
  
  
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
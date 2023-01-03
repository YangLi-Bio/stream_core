#' Calculate the distance for each enhancer to their nearest TSS,
#' which was borrowed from Signac
#'
#' @keywords internal
#'
DistanceToTSS <- function(peaks, genes, distance = 200000,
                          sep = c("-", "-")) {

  tss <- GenomicRanges::resize(x = genes, width = 1, fix = 'start') # find the TSS location
  genes.extended <- suppressWarnings(
    expr = Signac::Extend(
      x = tss, upstream = distance, downstream = distance
    )
  ) # extand the genomic range from the TSS till downstream/upstream 200000 bp
  overlaps <- GenomicAlignments::findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  ) # find the peaks overlaped with the extended genomic ranges of genes
  hit_matrix <- Matrix::sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  ) # build a sparse matrix to record the overlaps between peaks and extended genomic ranges of genes
  rownames(x = hit_matrix) <- Signac::GRangesToString(grange = peaks, sep = sep) # use peak names as the row names
  colnames(x = hit_matrix) <- genes.extended$gene_name # use gene names as the column names
  hit_matrix
}



#' Retain the longest transcripts for all exons
#'
#' @keywords internal
#' @import data.table
#'
CollapseToLongestTranscript <- function(ranges) {

  range.df <- data.table::as.data.table(x = ranges) # transform a GRanges object into a data frame
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
  gene.ranges <- GenomicRanges::makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = T # the information not used to build the data frame
    # will be retained in meta data
  )
  gene.ranges
}



#' Retain genes or enhancers which are within a nearby range of at least one enhancer or gene,
#' respectively
#'
#' @importFrom dplyr %>% filter
#'
#' @keywords internal
#'
filter_nearby_genes <- function(obj, distance = 500000, peak.assay = "ATAC") {

  # calculate the nearby genes
  gene.coords <- CollapseToLongestTranscript(ranges = Signac::Annotation(object = obj[[peak.assay]]))
  peaks <- Signac::StringToGRanges(rownames(obj[[peak.assay]])) # get peak coordinates
  distance.df <- Matrix::summary(DistanceToTSS(peaks = peaks, genes = gene.coords,
                                       distance = distance)) # distance matrix
  peak.names <- rownames(obj[[peak.assay]]) # get peak names
  gene.names <- gene.coords$gene_name # get gene names


  data.table::rbindlist(pbmcapply::pbmclapply(1:nrow(distance.df), function(i) {
    return(list(peak = peak.names[distance.df[i, 1]], gene = gene.names[distance.df[i, 2]]))
  }, mc.cores = min(parallel::detectCores(), nrow(distance.df))), fill = T) %>%
    dplyr::filter(gene %in% rownames(obj[["RNA"]]))
}



# LinksToGRanges <- function(linkmat, gene.coords, sep = c("-", "-")) {
#
#   tss <- GenomicRanges::resize(gene.coords, width = 1, fix = 'start')# the TSS
#   gene.idx <- sapply(
#     X = rownames(x = linkmat), # the gene names
#     FUN = function(x) {
#       which(x = x == tss$gene_name)[[1]]
#     }
#   )
#   tss <- tss[gene.idx] # select the TSS for each gene
#
#   # get midpoint of each peak
#   peak.ranges <- StringToGRanges(
#     regions = colnames(x = linkmat),
#     sep = sep
#   )
#   midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2) # the midpoint of the peaks
#
#   # convert to triplet form
#   dgtm <- as(object = linkmat, Class = "dgTMatrix")
#
#   # create data frame
#   df <- data.frame(
#     chromosome = as.character(x = seqnames(x = peak.ranges)[dgtm@j + 1]),
#     tss = start(x = tss)[dgtm@i + 1],
#     pk = midpoints[dgtm@j + 1],
#     score = dgtm@x,
#     gene = rownames(x = linkmat)[dgtm@i + 1],
#     peak = colnames(x = linkmat)[dgtm@j + 1]
#   )
#
#   # work out start and end coords (I did not understand these lines)
#   df$start <- ifelse(test = df$tss < df$pk, yes = df$tss, no = df$pk)
#   df$end <- ifelse(test = df$tss < df$pk, yes = df$pk, no = df$tss)
#   df$tss <- NULL
#   df$pk <- NULL
#
#   # convert to granges
#   gr.use <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
#   return(sort(x = gr.use))
# }




# link_signac <- function(x, distance = 500000,
#                         signac.score = 0,
#                         signac.pval = 1,
#                         min.cells = 10,
#                         peak.assay = 'ATAC') {
#
#   DefaultAssay(x) <- peak.assay # set 'ATAC' as he default assay
#   xx <- x # back up
#   x <- tryCatch(LinkPeaks(object = x, distance = distance,
#                            min.cells = min.cells,
#                            peak.assay = peak.assay, expression.assay = 'RNA',
#                            pvalue_cutoff = signac.pval,
#                            score_cutoff = signac.score, verbose = T),
#                 error = function(e) {
#                   0 }) # build linkages
#   if (!is.numeric(x)) {
#     signac.links <- data.frame(node1 = x[[peak.assay]]@links$peak,
#                                node2 = x[[peak.assay]]@links$gene,
#                                weight = x[[peak.assay]]@links$score)
#   } else {
#     x <- xx
#     signac.links <- filter_nearby_genes(obj = x) # link peaks to genes using heuristics
#     signac.links <- cbind(signac.links, pbmclapply(1 : nrow(signac.links),
#                                                    function(i) {
#       vx <- as.vector(x[["RNA"]][signac.links$gene[i]])
#       vy <- as.vector(x[[peak.assay]][signac.links$peak[i]])
#
#       if (sum(vx > 0) <= 0 | sum(vy > 0) <= 0) {
#         return(0)
#       } else if (sd(vx) == 0 | sd(vy) == 0) {
#         return(1)
#       } else {
#         return(cor(x = vx, y = vy, method = "pearson"))
#       }
#     }, mc.cores = min(detectCores(), nrow(signac.links))) %>% unlist)
#     colnames(signac.links) <- c("node1", "node2", "weight")
#     signac.links <- signac.links[signac.links$weight > signac.score,] # filter linkages
#   } # no peaks are linked to the genes
#   if (nrow(signac.links) < 1) {
#     return(NULL)
#   }
#
#
#   # Calculate weights
#   max.weight <- max(signac.links$weight)
#   min.weight <- min(signac.links$weight)
#   diff <- max.weight - min.weight
#   ifelse(diff > 0, signac.links$weight <- (max.weight - signac.links$weight) /
#            diff, signac.links$weight <- 0) # normalize the weights
#   message ('Identified ', nrow(signac.links), ' enhancer-gene linkages.\n')
#
#   signac.links
# }



# get_coherent_peak_gene_pairs <- function(peak_distance_matrix,
#                                          HBC.rna, HBC.atac) {
#
#   # Libraries
#   library(data.table)
#
#
#   row.col <- which(peak_distance_matrix > 0, arr.ind = T) # get the nonzero elements
#   peak.gene <- rbindlist(apply(row.col, 1, function(rr) {
#     return(list(rownames(peak_distance_matrix)[rr[1]], colnames(peak_distance_matrix)[rr[2]]))
#   })) # convert row and column ids into peaks and rows
#   colnames(peak.gene) <- c("peak", "gene")
#   peak.gene.weight <- cbind(peak.gene, weight = apply(peak.gene, 1, function(rr) {
#     return(length(which(HBC.atac[rr[1], ] > 0 & HBC.rna[rr[2], ] > 0)))
#   }))
#   peak.gene.weight <- peak.gene.weight[with(peak.gene.weight, order(gene, -weight))] # sort the data frame
#
#
#   return(tapply(peak.gene.weight$peak, peak.gene.weight$gene, function(i) {
#     head(i, n = 1) # select the pairs of peaks and genes
#   }))
# }



#' Load the genomic annotations of an organism
#'
#' @keywords internal
#'
load_database <- function(org = "hg38") {

  ifelse (grepl("^mm", org), enddb <- "org.Mm.eg.db",
          enddb <- "org.Hs.eg.db")


  return(enddb)
}



#' Convert ENSEMBL IDs to gene symbols
#'
#' @keywords internal
#'
ensembl_to_symbol <- function(ensembl.ll, org = org) {

  # library(AnnotationDbi)
  org.db <- load_database(org = org)
  # library(enddb, character.only = T)


  # Remove the version numbers
  key.ll <- ensembl.ll
  if (grepl("^mm", org)) {
    key.ll <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1", key.ll)
  } else {
    key.ll <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1", key.ll)
  }


  # Map Ensembl IDs to gene symbols
  require(org.db, character.only = T)
  symbol.ll <- AnnotationDbi::mapIds(x = get(org.db),
                      keys = key.ll,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")


  return(symbol.ll)
}



#' Load the annotation database of an organism
#'
#' @keywords internal
#'
load_annotation <- function(org = "hg38") {

  ifelse (grepl("^mm", org), enddb <- "EnsDb.Mmusculus.v75",
          enddb <- "EnsDb.Hsapiens.v75")


  enddb
}



#' Generate a \code{GRanges} object composed of gene annotations
#'
#' @keywords internal
#'
build_gene_GRanges <- function(org = "hg38") {

  enddb <- load_annotation(org = org)
  require(enddb, character.only = T) # if get the error of "cli", please run "update.packages("cli")"
  annotations <- Signac::GetGRangesFromEnsDb(get(enddb))
  ensembldb::seqlevelsStyle(annotations) <- 'UCSC'
  gene.gr <- CollapseToLongestTranscript(annotations)

  return(gene.gr)
}



#' Link enhancers to genes witrhin a distance cutoff
#'
#'@importFrom Matrix summary
#'
#' @keywords internal
#'
link_peaks_to_genes <- function(peak.obj = c("chrX-192989-220023", "chr2-178095031-178129859"),
                                gene.obj = c("PLCXD1", "NFE2L2"),
                                org = "hg38", distance = 250000) {

  peak.gr <- peak.obj
  gene.annotation <- build_gene_GRanges(org = org)
  if (grepl("^ENSG|ENSMUG", gene.annotation$gene_name[1])) {
    symbol.ll <- ensembl_to_symbol(ensembl.ll = gene.annotation$gene_name,
                                   org = org)
    gene.gr <- gene.annotation[which(gene.annotation$gene_name %in% symbol.ll)]
  } else {
    gene.gr <- gene.annotation[which(gene.annotation$gene_name %in% gene.obj)]
  }


  # Link peaks to genes
  if (is.numeric(distance)) {
    message ("Finding nearby genes for each peak within ", distance, " bp ...")
    summ <- DistanceToTSS(peaks = peak.gr, genes = gene.gr,
                          distance = distance, sep = c("-", "-")) %>% summary
  } else if (distance == "gene") { # to-do : find the closest gene for each peak using Signac
    message ("Finding the closest peak for each gene ...")
    summ <- Signac::nearest(gene.gr, peak.gr)
    summ <- summ[which(!is.na(summ))]
    summ <- data.frame(i = summ, j = seq_along(summ), x = rep(1, length(summ)))
  }
  cis.gr <- peak.gr[summ$i]
  GenomicRanges::mcols(cis.gr)$gene <- gene.gr$gene_name[summ$j]


  return(cis.gr)
}

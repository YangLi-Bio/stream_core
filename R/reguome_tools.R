# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
data.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/"


# # Source the R scripts
# source(paste0(code.dir, "TFBS_list.R"))


# Annotate CREs with TFs
#' @import Signac pbmcapply dplyr
find_TFBS <- function(m.atac, TFBS.list, org = "hg38") {

  # Libraries
  library(Signac)
  library(pbmcapply)
  library(dplyr)


  peaks <- rownames(m.atac)
  if (grepl("mm", org)) {
    jaspar.sites <- TFBS.list[["Mouse"]]
  } else {
    jaspar.sites <- TFBS.list[["Human"]]
  }
  overlap <- findOverlaps(query = StringToGRanges(peaks), subject = jaspar.sites$peak)
  if (length(overlap) < 1) {
    stop ("No TF is found to bind the CREs.\n")
  }
  TF.peak.df <- do.call(rbind.data.frame, pbmclapply(seq_along(overlap), function(x) {
    return(data.frame(TF = jaspar.sites$TF[overlap@to[x]], peak = peaks[overlap@from[x]]))
  }, mc.cores = min(detectCores(), length(overlap)))) %>% distinct # build the TF-peak data frame
  message ("Finished identified ", nrow(TF.peak.df), " TF-peak pairs.\n")
  peak.dfs <- split(x = TF.peak.df$TF, f = TF.peak.df$peak)
  # nest list, where the names are peaks and the elements are TF lists

  TF.dfs <- split(x = TF.peak.df$peak, f = TF.peak.df$TF)


  list(CRE = peak.dfs, TF = TF.dfs)
}


# Discover cell-subpopulation-active TF-target pairs
#' @import pbmcapply igraph
find_TF_gene <- function(G.list, bound.TFs,
             binding.CREs) {

  # Libraries
  library(pbmcapply)
  library(igraph)


  return(pbmclapply(G.list, function(G) {
    G.vs <- as_ids(V(G)) # all nodes
    peak.vs <- G.vs[grep("^chr", G.vs)] # peak nodes
    gene.vs <- setdiff(G.vs, peak.vs) # gene nodes
    g.ends <- as.data.frame(ends(graph = G, es = E(G)[peak.vs %--% gene.vs]))
    gene.peaks <- split(g.ends, g.ends$V2) %>% lapply(., `[[`, ("V1")) # gene-peak relations
    peak.genes <- split(g.ends, g.ends$V1) %>% lapply(., `[[`, ("V2")) # peak-gene relations
    TF.genes <- lapply(names(binding.CREs), function(t) {
      return(unique(unlist(peak.genes[intersect(binding.CREs[[t]], names(peak.genes))])))
    })
    names(TF.genes) <- names(binding.CREs)
    TF.genes <- TF.genes[!sapply(TF.genes, is.null)]


    # build gene.TF relations
    gene.TFs <- lapply(names(gene.peaks), function(g) {
      return(unique(unlist(bound.TFs[intersect(gene.peaks[[g]], names(bound.TFs))])))
    })
    names(gene.TFs) <- names(gene.peaks)
    gene.TFs <- gene.TFs[!sapply(gene.TFs, is.null)]


    return(list(TF.genes = TF.genes, gene.TFs = gene.TFs))
  }, mc.cores = detectCores()))
}

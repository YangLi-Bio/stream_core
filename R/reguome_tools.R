#' Annotate enhancers with TFs
#'
#' @importFrom dplyr %>% distinct
#'
#' @keywords internal
#'
find_TFBS <- function(m.atac, TFBS.list, org = "hg38") {

  peaks <- rownames(m.atac)
  if (grepl("mm", org)) {
    jaspar.sites <- TFBS.list[["Mouse"]]
  } else {
    jaspar.sites <- TFBS.list[["Human"]]
  }
  overlap <- GenomicAlignments::findOverlaps(query = Signac::StringToGRanges(peaks),
                                             subject = jaspar.sites$peak)
  if (length(overlap) < 1) {
    stop ("No TF is found to bind the enhancers.")
  }
  TF.peak.df <- do.call("rbind", pbmcapply::pbmclapply(seq_along(overlap), function(x) {
    return(data.frame(TF = jaspar.sites$TF[overlap@to[x]], peak = peaks[overlap@from[x]]))
  }, mc.cores = min(parallel::detectCores(), length(overlap)))) %>% dplyr::distinct()
  # build the TF-enhancer data frame

  message ("Finished identified ", nrow(TF.peak.df), " TF-enhancer pairs.")
  peak.dfs <- split(x = TF.peak.df$TF, f = TF.peak.df$peak)
  # nest list, where the names are peaks and the elements are TF lists

  TF.dfs <- split(x = TF.peak.df$peak, f = TF.peak.df$TF)


  list(CRE = peak.dfs, TF = TF.dfs)
}



#' Discover cell-subpopulation-active TF-target pairs
#'
#' @importFrom dplyr %>%
#' @importFrom igraph %--%
#'
#' @keywords internal
#'
find_TF_gene <- function(G.list, bound.TFs,
             binding.CREs) {

  return( pbmcapply::pbmclapply(G.list, function(G) {
    G.vs <- igraph::as_ids(V(G)) # all nodes
    peak.vs <- G.vs[grep("^chr", G.vs)] # peak nodes
    gene.vs <- setdiff(G.vs, peak.vs) # gene nodes
    g.ends <- as.data.frame(igraph::ends(graph = G, es = igraph::E(G)[peak.vs %--% gene.vs]))
    gene.peaks <- split(g.ends, g.ends$V2) %>% lapply(., `[[`, ("V1")) # gene-peak relations
    peak.genes <- split(g.ends, g.ends$V1) %>% lapply(., `[[`, ("V2")) # peak-gene relations
    TF.genes <- lapply(names(binding.CREs), function(t) {
      return(unique(unlist(peak.genes[intersect(binding.CREs[[t]],
                                                names(peak.genes))])))
    })
    names(TF.genes) <- names(binding.CREs)
    TF.genes <- TF.genes[!sapply(TF.genes, is.null)]


    # build gene.TF relations
    gene.TFs <- lapply(names(gene.peaks), function(g) {
      return(unique(unlist(bound.TFs[intersect(gene.peaks[[g]],
                                               names(bound.TFs))])))
    })
    names(gene.TFs) <- names(gene.peaks)
    gene.TFs <- gene.TFs[!sapply(gene.TFs, is.null)]


    return(list(TF.genes = TF.genes, gene.TFs = gene.TFs))
  }, mc.cores = parallel::detectCores()) )
}

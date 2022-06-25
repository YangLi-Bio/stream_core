# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
data.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/"


# Subset a Seurat object
#' @import Seurat pbapply
subset_object <- function(block.list, object, seed.ratio = 0, atac.assay = 'ATAC',
                          atac.dis, max.peaks = 3000, links.df, min.cells = 0) {

  # Libraries
  library(pbapply)
  library(Seurat)


  atac.dis <- atac.dis[intersect(unique(links.df[links.df$gene %in% rownames(object[["RNA"]])]$peak),
                                 rownames(atac.dis)),]

  return(pblapply(block.list, function(block) {
    n.cutoff <- 0
    r.sum <- rowSums(atac.dis[, block$cells]) # rowwise sums
    r.sum <- r.sum[which(r.sum > n.cutoff)] # remove the zero sums
    r.sum <- r.sum[order(r.sum, decreasing = T)] # rank the row sums
    block.peaks <- names(r.sum) # get the consistent peaks
    if (length(block.peaks) < 1) {
      return(NULL)
    }
    if (length(block.peaks) > max.peaks) {
      block.peaks <- block.peaks[1:max.peaks]
    }

    return(subset(x = object, features = c(rownames(object[['RNA']]), block.peaks),
                  cells = block$cells))
  }))
}


# Build heterogeneous graph
#' @import Signac cicero BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg38 BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm9 pbapply igraph
build_graph <- function(obj.list,
                        distance = 500000, cicero.covar = 0,
                        org.gs = BSgenome.Hsapiens.UCSC.hg38,
                        signac.score = 0, signac.pval = 0.5,
                        min.cells = 10,
                        peak.assay = 'ATAC') {

  # Libraries
  library(Signac)
  library(cicero)
  library(org.gs, character.only = T)
  library(pbapply)
  library(igraph)


  return(pblapply(obj.list, function(x) {
    signac.links <- link_signac(x, distance = distance,
                                signac.score = signac.score,
                                signac.pval = signac.pval,
                                min.cells = min.cells,
                                cell.weight = cell.weight,
                                peak.assay = peak.assay)
    if (is.null(signac.links)) {
      return(NULL)
    }


    # Link CREs to CREs
    x.atac <- GetAssayData(object = x, slot = 'data',
                           assay = atac.assay) # extract the ATAC assay
    ifelse (nrow(x.atac) <= 1000 | ncol(x.atac) <= 1000,
            cicero.links <- link_cor(x.atac, distance = distance,
                                     cicero.covar = cicero.covar),
            cicero.links <- tryCatch(link_cicero(x.atac,
                                                 distance = distance, org.gs = org.gs,
                                                 cicero.covar = cicero.covar),
                                     error = function(e) {
                                       message ('Error: something is wrong with Cicero.\n')
                                       0
                                     }))
    if (is.null(cicero.links) | "data.frame" %!in% class(cicero.links)) { # error exists
      cicero.links <- link_cor(x.atac, distance = distance,
                               cicero.covar = cicero.covar,
                               cell.weight = cell.weight)
    }


    quiet(ifelse ("data.frame" %in% class(cicero.links), links <- rbind(signac.links, cicero.links),
                  links <- signac.links)) # bind two data frames
    G <- graph_from_data_frame(links, directed = F) # build graph
    E(G)$weight <- links$weight # assign the edge weights


    return(G)
  }))
}


# Seeding based on Steiner forest problem (SFP) model
#' @import igraph pbmcapply
SFP_seeding <- function(block.list, G.list, obj.list, bound.TFs, binding.CREs,
                           TFGene.pairs, rna.dis, atac.dis, KL = 6, P = NULL,
                           Q = NULL, score.cutoff = 1) {

  # Libraries
  library(igraph)
  library(pbmcapply)


  seed.es <- Reduce(rbind, pbmclapply(seq_along(block.list), function(i) {
    G <- G.list[[i]] # the graph
    G.vs <- as_ids(V(G)) # all nodes
    peak.vs <- G.vs[grep("^chr", G.vs)] # peak nodes
    terminal <- intersect(block.list[[i]]$genes, G.vs) # gene nodes
    if (length(terminal) == 0) {
      return(NULL)
    }
    G <- induced_subgraph(graph = G, vids = c(peak.vs, terminal))
    es.df <- MST_to_SFP(G, terminal) %>% rearrange_cols # solve the SFP model using MST
    x.df <- remove_steiner(es.df) # remove the edges whose two nodes are peaks
    colnames(x.df) <- c('steiner_node', 'terminal_node') # name the columns
    x.df <- cbind(terminal = rep(i, nrow(x.df)), x.df)


    return(x.df)
  }, mc.cores = min(detectCores(), length(G.list))))
  seeds <- Reduce(c, pbmclapply(unique(seed.es$terminal), function(i) {
    G <- G.list[[i]]
    # G <- G.list
    block.cells <- block.list[[i]]$cells # the cells where genes are coexpressed
    df <- seed.es[seed.es$terminal == i, ] # get the data frame
    gene.TFs <- TFGene.pairs[[i]]$gene.TFs # gene-TF relations of this cell cluster
    TF.genes <- TFGene.pairs[[i]]$TF.genes # TF-gene relations of this cell cluster
    TFs <- Reduce(c, gene.TFs[unique(df$terminal_node)]) # get the regulating TFs
    names(TFs) <- NULL
    TF.freq <- table(TFs) # the occurring frequency of TFs
    TF.freq <- TF.freq[which(TF.freq > score.cutoff)] # remove the TFs regulating no gene
    if (length(TF.freq) < 2) {
      return(NULL)
    }

    TF.freq <- TF.freq[order(TF.freq, decreasing = T)] # filter TFs
    TF.top <- names(TF.freq[1:min(length(TF.freq), TOP_TFS)]) # select the top-ranked TFs
    ter.seeds <- lapply(TF.top, function(tt) {
      overlap.genes <- intersect(TF.genes[[tt]], unique(df$terminal_node)) # genes regulated by the TF
      if (length(overlap.genes) <= score.cutoff) {
        return(NULL)
      }

      overlap.peaks <- intersect(unique(df$steiner_node[df$terminal_node %in% overlap.genes]),
                                 binding.CREs[[tt]]) # peaks bound the the TF and meanwhile linked to
      if (length(overlap.peaks) < 1) {
        return(NULL)
      }
      overlap.genes <- df$terminal_node[df$terminal_node %in% overlap.genes &
                                          df$steiner_node %in% overlap.peaks]
      atac.ratio <- ifelse (length(overlap.peaks) > 1,
                            min(rowSums(atac.dis[overlap.peaks, block.cells, drop = F])) /
                              length(block.cells), # the cutoff of consistency between accessibility and expression
                            sum(atac.dis[overlap.peaks, block.cells]) / length(block.cells)
      )

      order.cells <- apply(rna.dis[overlap.genes, block.cells, drop = F], 2, sum) %>%
        sort(., decreasing = T) %>% names() # sort cells according to the sum of expression values
      HBC <- list(terminal = i, TF = tt, genes = overlap.genes, peaks = overlap.peaks,
                  cells = order.cells, atac.ratio = atac.ratio, score = 0, weight = 1)
      HBC$score <- score_HBC(HBC = HBC, KL = 6, m = rna.dis, P = P, Q = Q, G = G)


      return(HBC)
    })
  }, mc.cores = detectCores())) %>% discard(is.null)
  score.ll <- sapply(seeds, "[[", ("score")) # obtain the scores
  seeds <- seeds[order(score.ll, decreasing = T)] # sort the seeds in decreasing order of scores


  return(seeds[sapply(seeds, "[[", ("score")) > score.cutoff])
}

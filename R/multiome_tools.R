# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
data.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/"


# macros
VERY_SMALL <- 0.00001


# Subset a Seurat object
#' @import Seurat pbapply
subset_object <- function(block.list, object, seed.ratio = 0, peak.assay = 'ATAC',
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
  # library(org.gs)
  library(pbapply)
  library(igraph)


  return(pblapply(obj.list, function(x) {
    signac.links <- link_signac(x, distance = distance,
                                signac.score = signac.score,
                                signac.pval = signac.pval,
                                min.cells = min.cells,
                                cell.weight = cell.weight,
                                peak.assay = peak.assay)
    if (is.null(signac.links) | nrow(signac.links) < 1) {
      return(NULL)
    }


    # Link CREs to CREs
    x.atac <- GetAssayData(object = x, slot = 'data',
                           assay = peak.assay) # extract the ATAC assay
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
    if (is.null(cicero.links) | "data.frame" %!in%
        class(cicero.links)) { # error exists
      cicero.links <- link_cor(x.atac, distance = distance,
                               cicero.covar = cicero.covar,
                               cell.weight = cell.weight)
    }


    quiet(ifelse ("data.frame" %in% class(cicero.links),
                  links <- rbind(signac.links, cicero.links),
                  links <- signac.links)) # bind two data frames
    G <- graph_from_data_frame(links, directed = F) # build graph
    E(G)$weight <- links$weight # assign the edge weights


    return(G)
  }))
}


# Construct Steiner forest based on minimum spanning trees
#' @import igraph dplyr
MST_to_SFP <- function(G, ter) {

  # Libraries
  library(igraph)
  library(dplyr)


  comp.info <- components(graph = G, mode = "weak") # calculate the connected components
  splited.ter <- as.factor(comp.info$membership[ter]) %>% split(., f = factor(.)) # split terminals into sets
  Reduce(rbind, lapply(seq_along(splited.ter), function(i) {
    if (length(splited.ter[[i]]) < 2) { # only one gene
      v.edges <- incident(graph = G, v = names(splited.ter[[i]]), mode = "all") # get all neighbors
      return(as.data.frame(ends(graph = G, es = v.edges[which.min(v.edges$weight)])))
    } else if (length(splited.ter[[i]]) == 2) { # multiple genes
      return(as.data.frame(ends(graph = G,
                                es = shortest_paths(graph = G, from = names(splited.ter[[i]])[[1]],
                                                    to = names(splited.ter[[i]])[[2]], output = "epath",
                                                    weights = E(G)$weight) %>% .$epath %>%
                                  `[[` (1)))) # the shortest path
    } else {
      path <- shortest_paths(graph = G, from = names(splited.ter[[i]])[[1]],
                             to = names(splited.ter[[i]])[[2]], output = "both",
                             weights = E(G)$weight)
      path.vs <- as_ids(path$vpath[[1]]) # nodes
      for (j in 3 : length(splited.ter[[i]])) {
        path <- shortest_paths(graph = G, from = path.vs,
                               to = names(splited.ter[[i]])[[j]], output = "both",
                               weights = E(G)$weight)
        path.vs <- union(path.vs, as_ids(path$vpath[[1]])) # add the new nodes
      }
      sub.G <- induced_subgraph(graph = G, vids = path.vs) # build subgraph
      sub.mst <- mst(graph = sub.G,
                     weights = E(sub.G)$weight, algorithm = "prim") # get the minimum spanning tree


      return(as.data.frame(ends(graph = sub.mst, es = E(sub.mst))))
    }
  }))
}


# Remove steiner nodes, i.e., CREs
remove_steiner <- function(df) {
  Reduce(rbind, lapply(1:nrow(df), function(i) {
    if (grepl("^chr", df[i, 2])) {
      return(NULL)
    }


    return(df[i,])
  }))
}


# Score hybrid biclusters (HBCs)
#' @import dplyr
score_HBC <- function(HBC, m = NULL, KL = 6, Q = NULL,
                      P = NULL, G = NULL) {

  # Libraries
  library(dplyr)


  if (KL == 0) {
    R <- t(apply(m[HBC$genes, HBC$cells], 1, function(r) {
      n <- sum(r > 0)
      return(c(n, ncol(m[HBC$genes, HBC$cells]) - n))
    })) %>% `+` (VERY_SMALL) %>% `/` (ncol(m[HBC$genes, HBC$cells]) +
                                        2 * VERY_SMALL)
    C <- t(apply(m[HBC$genes, HBC$cells], 2, function(cc) {
      n <- sum(cc > 0)
      return(c(n, nrow(m[HBC$genes, HBC$cells]) - n))
    })) %>% `+` (VERY_SMALL) %>% `/` (nrow(m[HBC$genes, HBC$cells]) +
                                        2 * VERY_SMALL)


    return(mean(sapply(rownames(R), function(j) {
      return(R[j, 1] * log(R[j, 1] / Q[j, 1]) + R[j, 2] * log(R[j, 2] / Q[j, 2]))
    })) +
      mean(sapply(rownames(C), function(j) {
        return(C[j, 1] * log(C[j, 1] / P[j, 1]) + C[j, 2] * log(C[j, 2] / P[j, 2]))
      })))
  } else if (KL == 1) { # this score decrease gradually
    return(length(HBC$genes) * length(HBC$cells))
  } else if (KL == 2) { # this score decrease generally
    return(mean(E(G)[HBC$genes %--% HBC$peaks]$weight))
  } else if (KL == 3) { # this score decrease gradually
    return(length(E(G)[HBC$genes %--% HBC$peaks]) /
             length(HBC$genes) / length(HBC$peaks))
  } else if (KL == 4) { # this score decrease gradually
    return(sum(E(G)[HBC$genes %--% HBC$peaks]$weight) /
             length(HBC$genes) / length(HBC$peaks))
  } else if (KL == 5) { # this score decrease generally
    return(sum(E(G)[HBC$genes %--% HBC$peaks]$weight) /
             length(HBC$genes))
  }


  return(min(length(HBC$genes), length(HBC$cells)))
}


# Rearrange the columns of a data frame
rearrange_cols <- function(df) {

  Reduce(rbind, lapply(1 : nrow(df), function(i) {
    if (!grepl("^chr", df[i, 1])) {
      return(list(df[i, 2:1]))
    }

    return(df[i, ])
  }))
}


# Seeding based on Steiner forest problem (SFP) model
#' @import igraph pbmcapply dplyr scales
SFP_seeding <- function(block.list, G.list, obj.list, bound.TFs, binding.CREs,
                           TFGene.pairs, rna.dis, atac.dis, KL = 6, P = NULL,
                           Q = NULL, score.cutoff = 1) {

  # Libraries
  library(igraph)
  library(pbmcapply)
  library(dplyr)
  library(scales)


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
    block.cells <- block.list[[i]]$cells # the cells where genes are coexpressed
    df <- seed.es[seed.es$terminal == i,] # get the data frame
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
    TF.top <- names(TF.freq[1 : min(length(TF.freq), TOP_TFS)]) # select the top-ranked TFs
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
  }, mc.cores = detectCores()))
  seeds <- seeds[sapply(seeds, is.not.null)]
  score.ll <- sapply(seeds, "[[", ("score")) # obtain the scores
  seeds <- seeds[order(score.ll, decreasing = T)] # sort the seeds in decreasing order of scores


  return(seeds[sapply(seeds, "[[", ("score")) > score.cutoff])
}


# Get the list of RNA and ATAC matrices
#' @import Seurat
get_matrix_list <- function(m, obj.list, assay = "RNA") {

  # Libraries
  library(Seurat)


  return(lapply(obj.list, function(x) {
    return(m[intersect(rownames(x[[assay]]), rownames(m)),
             colnames(x[[assay]])])
  }))
}


# Whether the seed is eligible based on the similarity with other
# seeds related to the same TF
intra_eligible_seed <- function(s, h, gene.cutoff = 0.75) {

  cutoff <- gene.cutoff * min(length(s$genes), length(h$genes))
  overlap.genes <- intersect(s$genes, h$genes)
  if (length(overlap.genes) > cutoff) {
    return(F)
  }


  return(T)
}


# Get the top-ranked genes, 100 by default
#' @import dplyr igraph
get_top_genes_peaks <- function(HBC, cand.genes, cand.peaks,
                                top.ngenes = 100, G,
                                rna.m, atac.m) {

  # Libraries
  library(dplyr)
  library(igraph)


  GMULTI <- 1
  PMULTI <- 2 * GMULTI
  cand.genes <- setdiff(cand.genes, HBC$genes) # the genes to choose from
  if (length(cand.genes) < 1) {
    return(NULL)
  }
  cand.genes <- (rna.m[cand.genes, , drop = F] > 0) %>% rowSums %>%
    sort(., decreasing = T) %>%
    head(GMULTI * top.ngenes) %>% names # sort genes in decreasing
  cand.peaks <- setdiff(cand.peaks, HBC$peaks) %>% intersect(., rownames(atac.m))
  cand.peaks <- intersect(Reduce(union, lapply(cand.genes, function(v) {
    as_ids(neighbors(G, v, "all"))
  })), cand.peaks)
  if (length(cand.peaks) < 1) {
    return(NULL)
  }
  cand.peaks <- (atac.m[cand.peaks, , drop = F] > 0) %>% rowSums %>%
    sort(., decreasing = T) %>% head(PMULTI * top.ngenes) %>% names
  top.peaks.genes <- apply(ends(G, E(G)[cand.genes %--% cand.peaks]), 2, unique)

  if ("matrix" %in% class(top.peaks.genes)) {
    top.peaks.genes <- list(top.peaks.genes[, 1], top.peaks.genes[, 2])
  }


  return(top.peaks.genes)
}


# Expand the core part of the HBC
#' @import pbmcapply dplyr qualV
expand_core <- function(HBC, top.genes.peaks, rna.m, atac.m, G,
                        KL = F, Q = NULL, P = NULL, DELTA = 0.0015,
                        min.cells = 10) {

  # Libraries
  library(pbmcapply)
  library(dplyr)
  library(qualV)


  choice.genes <- top.genes.peaks[[2]] # genes to choose
  choice.peaks <- top.genes.peaks[[1]] # peaks to choose
  iter <- 0
  score.end <- F
  while (T) {
    if (length(choice.genes) < 1) {
      break
    }
    choice.info <- pbmclapply(choice.genes, function(g) {
      g.prof <- rna.m[g, HBC$cells] %>% sort(., decreasing = T)
      g.cells <- g.prof[which(g.prof > 0)] %>% names
      g.lcs <- LCS(g.cells, HBC$cells) # the longest common path
      g.cells <- g.lcs$LCS
      cell.cutoff <- g.lcs$LLCS * HBC$atac.ratio # the cutoff of cell numbers
      if (nrow(ends(G, E(G)[g %--% HBC$peaks])) > 0) {
        return(list(gene = g, peak = NA, cells = g.cells))
      }
      g.peaks <- as.data.frame(ends(G, E(G)[g %--% choice.peaks])) %>% `[[` (1) %>% unique

      if (length(g.peaks) < 1) {
        return(NULL)
      }
      p.cells <- atac.m[g.peaks, g.cells, drop = F] %>% rowSums # calculate the qualified peaks
      qual.peaks <- p.cells[which(p.cells >= cell.cutoff)] %>% names %>% head(n = 1) # select peaks

      if (length(qual.peaks) < 1) {
        return(NULL)
      }


      return(list(gene = g, peak = qual.peaks, cells = g.cells))
    }, mc.cores = detectCores()) %>% discard(is.null)
    qual.info <- choice.info[!is.na(choice.info %>% lapply(., "[[", ("gene")) %>% unlist)]
    if (length(qual.info) < 1) {
      break
    }
    add.id <- qual.info %>% lapply(., "[[", ("cells")) %>% sapply(., length) %>% which.max
    if (length(qual.info[[add.id]]$cells) < min.cells) {
      break
    }
    add.HBC <- HBC
    add.HBC$genes <- c(add.HBC$genes, qual.info[[add.id]]$gene)
    if (!is.na(qual.info[[add.id]]$peak)) {
      add.HBC$peaks <- c(add.HBC$peaks, qual.info[[add.id]]$peak)
    }
    add.HBC$cells <- qual.info[[add.id]]$cells
    add.HBC$score <- score_HBC(add.HBC, KL = 6, P = P, Q = Q, m = m, G = G)
    iter <- iter + 1
    if (add.HBC$score < HBC$score & length(HBC$genes) > 2) {
      break
    }
    choice.genes <- setdiff(choice.genes, qual.info[[add.id]]$gene)
    HBC <- add.HBC
  }


  return(HBC)
}


# Add other genes or peaks that are consistent with the HBC
#' @import dplyr igraph
expand_fuzzy <- function(HBC, G, cand.genes, cand.peaks, rna.m, atac.m,
                         g.cutoff = 1.0, m, mm,
                         P = NULL, Q = NULL, KL = F) {

  # Libraries
  library(dplyr)
  library(igraph)


  cand1.genes <- setdiff(cand.genes, HBC$genes) %>% intersect(., as_ids(V(G))) # remove the included genes
  cand1.peaks <- setdiff(cand.peaks, HBC$peaks) %>% intersect(., as_ids(V(G))) # remove the included peaks
  if (length(cand1.genes) < 1) {
    message(length(HBC$genes), "/", length(HBC$peaks), "/",
            length(HBC$cells), " genes/peaks/cells were included in this HBC.\n")

    return(HBC)
  }
  gene.cutoff <- g.cutoff * length(HBC$cells)
  cand2.genes <- cand1.genes[(rna.m[cand1.genes, HBC$cells, drop = F] > 0) %>%
                               rowSums >= gene.cutoff]
  if (length(cand2.genes) < 1) {
    message(length(HBC$genes), "/", length(HBC$peaks), "/",
            length(HBC$cells), " genes/peaks/cells were included in this HBC.\n")


    return(HBC)
  }
  retain.genes <- cand2.genes
  new.gene.links <- ends(G, E(G)[retain.genes %--% HBC$peaks]) %>%
    apply(., 2, unique) # new genes linked to old peaks
  if ("array" %in% class(new.gene.links)) { # in cases we got a matrix
    new.gene.links <- list(new.gene.links[, 1], new.gene.links[, 2])
  }
  if (length(new.gene.links) > 1) {
    add.genes <- c(HBC$genes, new.gene.links[[1]]) # add genes directly linked to old peaks
    retain.genes <- setdiff(retain.genes, new.gene.links[[1]]) # other genes to be added
  }
  extreme.cutoff <- HBC$atac.ratio
  peak.cutoff <- length(HBC$cells) * extreme.cutoff
  retain.peaks <- cand1.peaks[atac.m[cand1.peaks, HBC$cells, drop = F] %>% rowSums >= peak.cutoff]
  cand.es <- E(G)[retain.genes %--% retain.peaks]
  if (length(cand.es) < 1) {
    return(HBC)
  }
  cand.es.df <- cbind(as.data.frame(ends(graph = G, es = cand.es)), cand.es$weight) # convert to a data frame
  colnames(cand.es.df) <- c("peak", "gene", "weight")
  add.genes <- c(HBC$genes, unique(cand.es.df$gene)) # add new genes
  add.peaks <- c(HBC$peaks, unique(cand.es.df$peak[tapply(seq_along(cand.es.df$weight),
                                                          cand.es.df$gene, min)]))
  cc.g.cutoff <- g.cutoff * length(HBC$genes)
  cc.p.cutoff <- length(HBC$peaks) * extreme.cutoff
  cand.cells <- setdiff(colnames(m), HBC$cells) # non-included cells
  test.cells <- cand.cells[which(colSums(m[HBC$genes, cand.cells, drop = F]) >= cc.g.cutoff)]
  pri.add.cells <- test.cells[which(colSums(mm[HBC$peaks, test.cells, drop = F]) >= cc.p.cutoff)]
  add.cells <- c(HBC$cells, pri.add.cells)
  message(length(add.genes), "/", length(add.peaks), "/",
          length(add.cells), " genes/peaks/cells were included in this HBC.\n")
  HBC$genes <- add.genes
  HBC$peaks <- add.peaks
  HBC$cells <- add.cells
  HBC$score <- score_HBC(HBC, KL = 6, Q = Q, P = P, m = m, G = G) # score the HBC


  return(HBC)
}


# Expand a hybrid bicluster (HBC)
expand_HBC <- function(HBC, cand.genes, cand.peaks,
                       rna.m, atac.m, dual = F, G,
                       top.ngenes = 100, c.cutoff = 1.0,
                       KL = F, m = NULL, mm = NULL,
                       Q = NULL, P = NULL, min.cells = 10) {

  top.genes.peaks <- get_top_genes_peaks(HBC = HBC, cand.peaks = cand.peaks,
                                         cand.genes = cand.genes,
                                         G = G, rna.m = rna.m, atac.m = atac.m,
                                         top.ngenes = top.ngenes)
  if (is.null(top.genes.peaks)) {
    return(HBC)
  }
  if ("character" %in% class(top.genes.peaks)) {
    tmp.ll <- top.genes.peaks
    top.genes.peaks <- c()
    top.genes.peaks[[1]] <- tmp.ll[[1]]
    top.genes.peaks[[2]] <- tmp.ll[[2]]
  }
  HBC <- expand_core(HBC = HBC, top.genes.peaks = top.genes.peaks, rna.m = rna.m, atac.m = atac.m,
                     G = G, KL = KL, Q = Q, P = P, min.cells = min.cells) # expand the core part of the HBC
  HBC <- expand_fuzzy(HBC = HBC, G = G, cand.genes = cand.genes, cand.peaks = cand.peaks,
                      m = m, mm = mm, rna.m = rna.m, atac.m = atac.m, g.cutoff = c.cutoff, KL = KL,
                      P = P, Q = Q)


  return(HBC)
}


# Check whether a seed exists
exist_seed <- function(s, HBCs, same.terminal = T,
                       intra.cutoff = 0.50, inter.cutoff = 0.25,
                       peak.cutoff = 0.25) {

  for (h in HBCs) {

    # Filter seeds based on terminal
    if (same.terminal & s$terminal == h$terminal) {
      return(T)
    }


    # Filter seeds based on similarity
    if (s$TF == h$TF) {
      if (!intra_eligible_seed(s = s, h = h,
                               gene.cutoff = intra.cutoff)) {
        return(F)
      }
    } else {
      if (!inter_eligible_seed(s = s, h = h,
                               gene.cutoff = inter.cutoff,
                               peak.cutoff = peak.cutoff)) {
        return(F)
      }
    }
  }


  return(T)
}


# Check whether a seed is eligble based on the similarity with other
# seeds related to different TFs
inter_eligible_seed <- function(s, h,
                                gene.cutoff = 0.25,
                                peak.cutoff = 0.25) {

  # check overlaped genes
  cutoff <- gene.cutoff * min(length(s$genes), length(h$genes)) # the cutoff of gene number
  overlap.genes <- intersect(s$genes, h$genes)
  if (length(overlap.genes) <= cutoff) {
    return(T)
  }
  cutoff <- peak.cutoff * min(length(s$peaks), length(h$peaks)) # the cutoff of gene number
  overlap.peaks <- intersect(s$peaks, h$peaks)
  if (length(overlap.peaks) <= cutoff) {
    return(T)
  }


  return(F)
}


# Hybrid biclustering
#' @import igraph BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg38 BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm9 dplyr
hybrid_biclust <- function(seeds = seeds, rna.list = rna.list, atac.list = atac.list,
               top.ngenes = top.ngenes, bound.TFs = bound.TFs,
               binding.CREs = binding.CREs, G.list = G.list, TFGene.pairs = TFGene.pairs,
               c.cutoff = c.cutoff, KL = KL, org.gs = org.gs,
               rna.dis = rna.dis, atac.dis = atac.dis, min.cells = min.cells) {

  # Libraries
  library(igraph)
  # library(org.gs, character.only = T)
  library(dplyr)


  # Hybrid biclustering
  HBCs <- list() # the set of HBCs in nested list, which is the same as that of seeds
  id <- 0
  idd <- 0
  Q <- NULL
  P <- NULL
  if (KL == 0) {

    # proportion of nonzero and zero values in each row
    Q <- t(apply(rna.dis, 1, function(r) {
      n <- sum(r > 0)
      return(c(n, ncol(rna.dis) - n))
    })) %>% `+` (VERY_SMALL) %>% `/` (ncol(rna.dis) + 2 * VERY_SMALL)


    # proportion of nonzero and zero values in each column
    P <- t(apply(rna.dis, 2, function(cc) {
      n <- sum(cc > 0)
      return(c(n, nrow(rna.dis) - n))
    })) %>% `+` (VERY_SMALL) %>% `/` (nrow(rna.dis) + 2 * VERY_SMALL)
  }


  for (HBC in seeds) {
    idd <- idd + 1
    message ("Processing the seed : ", idd, " ...\n")
    if (!exist_seed(s = HBC, HBCs = HBCs, same.terminal = same.terminal,
                    intra.cutoff = intra.cutoff,
                    inter.cutoff = inter.cutoff,
                    peak.cutoff = peak.cutoff)) {
      next
    }
    id <- id + 1


    message("Processing the ", id, " hybrid bicluster ...\n")
    HBC <- expand_HBC(HBC = HBC, cand.peaks = binding.CREs[[HBC$TF]],
                      rna.m = rna.list[[HBC$terminal]], atac.m = atac.list[[HBC$terminal]],
                      dual = dual, cand.genes = TFGene.pairs[[HBC$terminal]]$TF.genes[[HBC$TF]],
                      G = G.list[[HBC$terminal]], m = rna.dis, mm = atac.dis,
                      top.ngenes = top.ngenes, c.cutoff = c.cutoff,
                      KL = KL, Q = Q, P = P, min.cells = min.cells)
    HBC$seed <- idd
    HBCs[[id]] <- HBC
  }


  return(HBCs)
}


# Merge significantly overlapped HBCs
#' @import dplyr pbmcapply data.table
merge_HBCs <- function(HBCs, stat = T, phyper.cutoff = 0.05,
                       rna.dis, atac.dis) {

  # Libraries
  library(dplyr)
  librry(pbmcapply)
  library(data.table)


  if (length(HBCs) < 2) {
    return(HBCs)
  }
  HBCs <- HBCs[HBCs %>% sapply(., "[[", ("genes")) %>% sapply(., length) > 2]
  HBCs <- HBCs[HBCs %>% sapply(., "[[", ("cells")) %>% sapply(., length) > 2]
  sorted.HBCs <- HBCs[order(HBCs %>% sapply(., "[[", ("score")), decreasing = T)]
  HBC.TFs <- unique(sapply(sorted.HBCs, "[[", "TF"))
  new.HBCs <- do.call("c", pbmclapply(HBC.TFs, function(tf) {
    TF.HBCs <- sorted.HBCs[sapply(sorted.HBCs, "[[", "TF") == tf]
    if (length(TF.HBCs) < 2) {
      return(TF.HBCs)
    }
    if (stat) {
      TF.cells <- Reduce(union, unlist(sapply(TF.HBCs, "[[", "cells")))
      N <- length(TF.cells)
      comb.pairs <- combn(seq_along(TF.HBCs), 2)
      adj.p.cutoff <- phyper.cutoff * length(TF.HBCs) * length(TF.HBCs)
      tri.df <- rbindlist(lapply(1:ncol(comb.pairs), function(j) {
        l1 <- comb.pairs[1, j]
        l2 <- comb.pairs[2, j]
        h1 <- TF.HBCs[[l1]]$cells
        h2 <- TF.HBCs[[l2]]$cells
        m <- length(h1)
        k <- length(h2)
        q <- length(intersect(h1, h2))
        weight <- phyper(q - 1, m, N - m, k, lower.tail = F)
        return(list(node1 = l1, node2 = l2, weight = weight))
      }), fill = T) %>% dplyr::filter(weight <= adj.p.cutoff) %>%
        dplyr::select(node1, node2) # build an igraph object
      if (nrow(tri.df) < 1) { # no need to merge
        return(TF.HBCs)
      }
      TF.cliques <- max_cliques(graph_from_data_frame(d = tri.df, directed = F)) # find maximal cliques
      merged.TF.HBCs <- lapply(TF.cliques, function(k) {
        HBC.list <- as.numeric(as_ids(k)) # the terminals to obtain the merged HBC
        clique.HBCs <- TF.HBCs[HBC.list] # the HBCs belonging to this clique
        terminal <- as.vector(sapply(clique.HBCs, "[[", "terminal"))
        genes <- unique(unlist(lapply(clique.HBCs, "[[", "genes")))
        peaks <- unique(unlist(lapply(clique.HBCs, "[[", "peaks")))
        cells <- unique(unlist(lapply(clique.HBCs, "[[", "cells")))
        atac.ratio <- mean(apply(atac.dis[peaks, cells], 1, sum)) / length(cells)
        new.HBC <- list(terminal = terminal, TF = tf, genes = genes, peaks = peaks,
                        cells = cells, atac.ratio = atac.ratio, score = 0)
        new.HBC$score <- score_HBC(new.HBC, KL = 6)


        new.HBC
      })
      ifelse(length(TF.cliques) < 2, covered.HBCs <- as.numeric(as_ids(TF.cliques[[1]])),
             covered.HBCs <- unique(unlist(sapply(TF.cliques, function(y) {
               as.numeric(as_ids(y))
             })))) # the set of merged HBCs
      isolated.TF.HBCs <- TF.HBCs[setdiff(seq_along(TF.HBCs), covered.HBCs)] # the isolated HBC
      return(c(merged.TF.HBCs, isolated.TF.HBCs)) # merge the new HBCs
    }
  }, mc.cores = detectCores()))


  return(new.HBCs[order(sapply(new.HBCs, "[[", "score"), decreasing = T)])
}


# Patch each HBC to add more genes and peaks
#' @import Seurat pbmcapply dplyr Signac
patch_HBCs <- function(merged.HBCs, binding.CREs, x, peak.ratio = NULL,
                       peak.assay = "ATAC", distance = 500000) {

  # Libraries
  library(Seurat)
  library(pbmcapply)
  library(dplyr)
  library(Signac)


  rna.m <- binarize(GetAssayData(x, slot = "data", assay = "RNA"))
  gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = x[[peak.assay]]))
  atac.m <- binarize(GetAssayData(x, slot = "data", assay = peak.assay))
  patched.HBCs <- pbmclapply(seq_along(merged.HBCs), function(i) {
    HBC <- merged.HBCs[[i]]
    HBC.rna <- rna.m[setdiff(rownames(rna.m), HBC$genes), HBC$cells, drop = F]
    genes.keep <- names(which(rowSums(HBC.rna) == ncol(HBC.rna))) %>%
      intersect(gene.coords$gene_name) # genes to keep
    if (length(genes.keep) < 1) { # No gene has coordinate annotation
      add.HBC <- HBC
    } else if (!is.finite(distance)) { # Without constraints of peaks
      HBC.seqs <- unique(seqnames(gene.coords[gene.coords$gene_name %in% HBC$genes]))
      genes.keep.coords <- gene.coords[gene.coords$gene_name %in% genes.keep] # coordinates of genes to add
      add.HBC <- HBC
      add.HBC$genes <- c(add.HBC$genes, genes.keep.coords[seqnames(genes.keep.coords) %in% HBC.seqs]$gene_name)
      add.HBC$score <- score_HBC(HBC = add.HBC, KL = 6)
    } else {
      HBC.gene.coords <- gene.coords[which(gene.coords$gene_name %in% genes.keep)] # coordinates to keep
      peak.HBC.distance <- DistanceToTSS(peaks = StringToGRanges(HBC$peaks), genes = HBC.gene.coords,
                                         distance = distance) # get the distance between genes and included peaks
      add.genes <- c(names(which(colSums(peak.HBC.distance) > 0))) # directly add the genes
      if (length(add.genes) < 1) {
        add.HBC <- HBC
      } else {
        HBC.gene.coords <- HBC.gene.coords[HBC.gene.coords$gene_name %!in% add.genes]
        HBC.atac <- atac.m[setdiff(binding.CREs[[HBC$TF]], HBC$peaks), HBC$cells, drop = F]
        peak.ratio <- quantile(rowSums(atac.m[HBC$peaks, HBC$cells,
                                              drop = F]))[[3]] / ncol(HBC.atac) # 25% quantile
        peaks.keep <- names(which(rowSums(HBC.atac) >= peak.ratio * ncol(HBC.atac))) # peaks to keep
        if (length(peaks.keep) < 1) { # Only add the genes linked to the peaks incorporated in the HBC
          add.HBC <- HBC
          add.HBC$genes <- c(add.HBC$genes, add.genes)
          add.HBC$score <- score_HBC(HBC = add.HBC, KL = 6)
        } else { # Have qualified peaks
          peaks.keep.GR <- StringToGRanges(peaks.keep) # get GRange objects
          peak_distance_matrix <- DistanceToTSS(peaks = peaks.keep.GR, genes = HBC.gene.coords,
                                                distance = distance) # get the distance between peaks and genes
          if (sum(peak_distance_matrix) == 0) {
            add.HBC <- HBC
            add.HBC$genes <- c(add.HBC$genes, add.genes)
            add.HBC$score <- score_HBC(HBC = add.HBC, KL = 6)
          } else {
            colnames(peak_distance_matrix) <- HBC.gene.coords$gene_name # rename the columns
            row.col <- get_coherent_peak_gene_pairs(peak_distance_matrix = peak_distance_matrix,
                                                    HBC.rna = HBC.rna, HBC.atac = HBC.atac)
            add.HBC <- HBC
            add.HBC$genes <- c(HBC$genes, names(row.col)) # update the genes
            add.HBC$peaks <- union(HBC$peaks, row.col) # update the peaks
            add.HBC$score <- score_HBC(HBC = add.HBC, KL = 6)
          }
        }
      }
    }
    add.gene.cells <- names(which(colSums(rna.m[HBC$genes,
                                                setdiff(colnames(rna.m),
                                                        HBC$cells),
                                                drop = F]) >= length(HBC$genes)))
    if (is.finite(distance)) {
      peak.ratio <- quantile(colSums(atac.m[HBC$peaks, HBC$cells,
                                            drop = F]))[[3]] / length(HBC$peaks)
    } else {
      peak.ratio <- quantile(colSums(atac.m[HBC$peaks, HBC$cells,
                                            drop = F]))[[4]] / length(HBC$peaks)
    }
    peak.cell.cutoff <- peak.ratio * length(HBC$peaks)
    add.peak.cells <- names(which(colSums(atac.m[HBC$peaks, add.gene.cells,
                                                 drop = F]) >= peak.cell.cutoff))
    add.HBC$cells <- c(add.HBC$cells, add.peak.cells)
    add.HBC$score <- score_HBC(add.HBC, KL = 6) # score the HBC


    return(add.HBC)
  }, mc.cores = min(detectCores(), length(merged.HBCs)))


  return(patched.HBCs[!sapply(patched.HBCs, is.null)])
}


# Calculate the pairwise similarity between eGRNs
#' @import pbmcapply dplyr
compute_original_sim <- function(HBCs, features = "genes") {

  # Libraries
  library(pbmcapply)
  library(dplyr)


  return(Reduce(rbind, pbmclapply(seq_along(HBCs), function(i) {
    sq1 <- length(HBCs[[i]][[features]]) * length(HBCs[[i]]$cells)

    sq1.sim <- lapply(seq_along(HBCs), function(j) {
      if(i < j) {
        return(NA)
      } else if (i == j) {
        return(list(HBC1 = i, HBC2 = i, Sim = 1, Same.TF = T))
      }

      sq2 <- length(HBCs[[j]][[features]]) * length(HBCs[[j]]$cells)
      sq <- length(intersect(HBCs[[i]][[features]], HBCs[[j]][[features]])) *
        length(intersect(HBCs[[i]]$cells, HBCs[[j]]$cells))
      sim <- sq / min(sq1, sq2) # original similarity

      ifelse(HBCs[[i]]$TF == HBCs[[j]]$TF & HBCs[[i]]$TF != 'NA',
             return(list(HBC1 = i, HBC2 = j, Sim = sim,
                         Same.TF = T)),
             return(list(HBC1 = i, HBC2 = j, Sim = sim,
                         Same.TF = F)))
    })

    return(rbindlist(sq1.sim[!is.na(sq1.sim)], fill = T))
  }, mc.cores = min(detectCores(), length(HBCs)))) %>%
    dplyr::filter(HBC1 != HBC2))
}


# Normalize the pairwise similarity
normalize_sim <- function(sim.df, HBCs) {

  same.df <- sim.df[sim.df$Same.TF, ] # pairs regulated by the same TF
  diff.df <- sim.df[!sim.df$Same.TF, ] # pairs regulated by different TFs
  norm.same.df <- same.df
  norm.same.df$Sim <- (norm.same.df$Sim - mean(norm.same.df$Sim)) / sd(norm.same.df$Sim)
  norm.diff.df <- diff.df
  norm.diff.df$Sim <- (norm.diff.df$Sim - mean(norm.diff.df$Sim)) / sd(norm.diff.df$Sim)
  return(rbind(norm.same.df, norm.diff.df, data.frame(HBC1 = seq_along(HBCs),
                                                      HBC2 = seq_along(HBCs),
                                                      Sim = rep(1, length(HBCs)),
                                                      Same.TF = rep(T, length(HBCs)))))

}


# Calculate pairwise similarity between HBCs
#' @import Matrix
compute_sim <- function(HBCs) {

  # Libraries
  library(Matrix)


  message ("Calculating the pairwise similarity between HBCs ...\n")
  gene.sim.df <- compute_original_sim(HBCs = HBCs) # compute the original similarity overlaps on expression level
  peak.sim.df <- compute_original_sim(HBCs = HBCs, features = "peaks")
  norm.gene.df <- normalize_sim(sim.df = gene.sim.df, HBCs = HBCs) # normalize the similarity matrix
  norm.peak.df <- normalize_sim(sim.df = peak.sim.df, HBCs = HBCs) # normalize the similarity matrix
  norm.sim.df <- norm.gene.df
  norm.sim.df$Sim <- sapply(1:nrow(norm.gene.df), function(i) {
    min(norm.gene.df$Sim[i], norm.peak.df$Sim[i])
  })
  sim.m <- sparseMatrix(i = norm.sim.df$HBC1, j = norm.sim.df$HBC2,
                        x = norm.sim.df$Sim) # generate a similarity matrix


  sim.m
}


# Calculate the submodular optimization value after adding a HBC
compute_add <- function(query, hit, sim.m, HBC.max.sim) {

  return(unlist(sapply(query, function(i) {
    add.sim <- max(sim.m[i, tail(hit, n = 1)],
                   sim.m[tail(hit, n = 1), i])
    ifelse (HBC.max.sim[[i]] < add.sim, return(add.sim),
            return(HBC.max.sim[[i]]))
  })))
}


# Select the next HBC
#' @import pbmcapply
select_HBC <- function(HBCs, sim.m, HBC.flag, HBC.max.sim,
                       regulon.ids) {

  # Libraries
  library(pbmcapply)


  max.sim.ll <- pbmclapply(seq_along(HBCs), function(x) {
    if (!HBC.flag[x]) {
      return(-1)
    } # do not consider the HBCs that have been already selected
    return(compute_add(query = seq_along(HBCs), hit = c(regulon.ids, x),
                       sim.m = sim.m,
                       HBC.max.sim = HBC.max.sim))
    # compute the objective value
  }, mc.cores = detectCores())
  which.HBC <- which.max(unlist(sapply(max.sim.ll, sum)))


  return(list(which.HBC, max.sim.ll[[which.HBC]]))
}


# Calculate the optimization function value
#' @import pbmcapply
calculate_sil <- function(selected.HBCs, links.df) {

  # Libraries
  library(pbmcapply)


  use.genes <- Reduce(union, lapply(selected.HBCs, "[[", "genes")) # all genes included in all regulons
  use.peaks <- Reduce(union, lapply(selected.HBCs, "[[", "peaks")) # all peaks included in all regulons
  use.links <- links.df[links.df$gene %in% use.genes &
                          links.df$peak %in% use.peaks, , drop = F] # useful links
  inner.links <- sum(unlist(pbmclapply(1:nrow(use.links), function(i) {
    whether.include <- F
    for (j in seq_along(selected.HBCs)) {
      if (use.links[i, 1] %in% selected.HBCs[[j]]$peaks &
          use.links[i, 2] %in% selected.HBCs[[j]]$genes) {
        whether.include <- T # include
        break
      }
    }
    return(whether.include)
  }, mc.cores = detectCores())))
  outer.links <- nrow(use.links) - inner.links # links not included
  sil.score <- inner.links - outer.links # Silhouette score
  message ("eGRN Silhouette score: ", sil.score, "\n.")


  return(sil.score) # silhouette score
}


# Perform submodular optimization
sub_mod <- function(HBCs, sim.m, G.list, n.cells, rna.list, block.list,
                    obj, peak.assay = "ATAC",
                    min.eGRNs = 100, submod.mode = 1,
                    distance = 500000) {

  regulon.ids <- c() # the list of HBC ids
  obj.list <- c() # the list of objective values
  HBC.flag <- rep(T, length(HBCs))
  HBC.max.sim <- rep(-1, length(HBCs))
  max.n <- 0 # the number of HBCs to select which yields the maximum objective value
  max.obj <- -2 # the current maximum objective value
  links.df <- filter_nearby_genes(obj, distance = distance,
                                  peak.assay = peak.assay)
  while (1) {
    double.output <- select_HBC(HBCs = HBCs, sim.m = sim.m,
                                HBC.flag = HBC.flag,
                                HBC.max.sim = HBC.max.sim,
                                regulon.ids = regulon.ids) # select the next HBC
    add.id <- double.output[[1]]
    HBC.max.sim <- double.output[[2]]
    regulon.ids <- c(regulon.ids, add.id) # add the new HBC
    message ("Finished evaluating the ", length(regulon.ids), " HBCs.\n")
    HBC.flag[add.id] <- F # mask this selected HBC
    obj.list <- c(obj.list, calculate_sil(selected.HBCs = HBCs[regulon.ids],
                                          links.df = links.df))
    if (length(regulon.ids) >= length(HBCs)) {
      break
    }
  }

  max.n <- tail(which(obj.list == max(obj.list[(min.eGRNs + 1) : length(HBCs)])), n = 1)
  message (max.n, " regulons yield the maximum score.\n")
  regulons <- HBCs[regulon.ids[1:max.n]] # select regulons


  return(list(regulons = regulons, obj = obj.list))
}


# Add putative genes to an eGRN if required
#' @import Seurat Signac pbapply
expand_eGRNs <- function(obj, submod.HBCs, peak.assay = 'ATAC',
                         distance = Inf) {

  # Libraries
  library(Seurat)
  library(Signac)
  librarY(pbapply)


  # Expansion
  rna.m <- GetAssayData(object = obj, slot = "data", assay = "RNA") > 0 # get the RNA expression matrix
  atac.m <- GetAssayData(object = obj, slot = "data", assay = peak.assay) > 0
  gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = obj[[peak.assay]])) # coordinates
  distance.df <- DistanceToTSS(peaks = StringToGRanges(rownames(atac.m)),
                               genes = gene.coords,
                               distance = distance) # distance matrix
  expanded.eGRNs <- pblapply(submod.HBCs, function(x) {
    add.genes <- setdiff(names(which(apply(rna.m[, x$cells], 1, sum) >=
                                       length(x$cells) * expand.cutoff)),
                         x$genes) # select additional genes
    add.coords <- gene.coords[gene.coords$gene_name %in% add.genes] # subset the genes to add
    overlap.genes <- intersect(add.genes, colnames(distance.df)) # some genes have not coordinates
    add.dist <- summary(distance.df[, overlap.genes]) # get the distance matrix
    x$gene.status <- c(rep(T, length(x$genes)), rep(F, length(add.genes))) # record the status of genes
    x$genes <- c(x$genes, add.genes) # extend genes
    x
  })


  expanded.eGRNs
}

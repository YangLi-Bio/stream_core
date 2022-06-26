# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
data.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/"


# macros
VERY_SMALL <- 0.00001


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


# Seeding based on Steiner forest problem (SFP) model
#' @import igraph pbmcapply dplyr
SFP_seeding <- function(block.list, G.list, obj.list, bound.TFs, binding.CREs,
                           TFGene.pairs, rna.dis, atac.dis, KL = 6, P = NULL,
                           Q = NULL, score.cutoff = 1) {

  # Libraries
  library(igraph)
  library(pbmcapply)
  library(dplyr)


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
  }, mc.cores = detectCores())) %>% discard(is.null)
  score.ll <- sapply(seeds, "[[", ("score")) # obtain the scores
  seeds <- seeds[order(score.ll, decreasing = T)] # sort the seeds in decreasing order of scores


  return(seeds[sapply(seeds, "[[", ("score")) > score.cutoff])
}


# Get the list of RNA and ATAC matrices
#' @import Seurat
get_matrix_list <- function(m, obj.list, assay = "RNA") {

  # Libraries
  library(Seurat)


  return(lapply(obj.list, function(obj) {
    return(m[rownames(obj[[assay]]), colnames(obj[[assay]])])
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
  library(org.gs, character.only = T)
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

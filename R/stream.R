#' @include cistrome_tools.R
#' @include epigenome_tools.R
#' @include multiome_tools.R
#' @include reguome_tools.R
#' @include transcriptome_tools.R
#' @include utilities.R
#'
NULL

#' Identify enhancer regulons (eRegulons) from jointly profiled scRNA-seq and scATAC-seq data
#'
#' @param obj A \code{Seurat} object composed of both scRNA-sed and scATAC-seq assays
#' @param top.peaks The number of top-ranked peaks to identify the core part of hybrid biclusters (HBCs),
#' 3000 by default
#' @param min.cells The cutoff of minimum number of cells for quality control (QC), 10 by default
#' @param out.dir The directory to save the intermediate or final results, "./" by default
#' @param org The organism version, hg38 by default
#' @param top.ngenes The number of genes composing the core part of an HBC, 5 by default
#' @param c.cutoff The cutoff of consistency during hybrid biclustering process, 1.0 by default
#' @param n.blocks The cutoff of the maximum number of blocks output by IRIS-FGM, 100 by default
#' @param min.eRegs The minimum number of enhancer regulons (eRegulons) to output, 100 by default
#' @param peak.assay The name of the scATAC-seq assay, "ATAC" by default
#' @param distance The distance cutoff to build enhancer-enhancer relations, 250000 by default
#' @param BlockOverlap The cutoff of maximum overlap between blocks output by \code{IRIS-FGM}, 0.500 by default
#' @param Extension The consistency level to expand a block by \code{IRIS-FGM}, 0.70 by default
#' @param intra.cutoff The cutoff to calculate pairwise similarity among HBCs associated with the same TFs,
#' 1.0 by default
#' @param inter.cutoff The cutoff to compute pairwise similarity among genes in HBCs associated with different TFs,
#' 0.80 by default
#' @param peak.cutoff The cutoff to quantify pairwise similarity among enhancers in HBCs, 0.80 by default
#' @param var.genes The number of highly variable genes to predict used to identify the core part of HBCs,
#' 3000 by default
#' @param KL Which method to use for measuring the score of HBCs, "min.exp" by default, i.e.,
#' the smaller one between the numbers of genes and cells in a HBC
#' @param quantile.cutoff The quantile cutoff of the ratio of HBC cells, where enhancers are accessible,
#' 4 by default, indicating the top-25% among ranks
#' @param submod.step The step size of the number of HBCs for submodular optimization, 30 by default
#'
#' @rdname run_stream
#' @concept run_stream
#'
#' @importFrom dplyr %>%
#'
#'
#' @export
#' @return When running on a \code{Seurat} object,
#' returns a list of enhancer regulons, i.e., eRegulons, saved in a nested list
run_stream <- function(obj = NULL,
                       peak.assay = "ATAC",
                       var.genes = 3000,
                       top.peaks = 3000,
                       min.cells = 10,
                       out.dir = "./",
                       org = "hg38",
                       top.ngenes = 5,
                       c.cutoff = 1.0,
                       n.blocks = 100,
                       distance = 250000,
                       ifWeighted = TRUE,
                       cicero.covar = -Inf,
                       signac.score = -Inf,
                       signac.pval = Inf,
                       intra.cutoff = 1.0,
                       inter.cutoff = 0.80,
                       peak.cutoff = 0.80,
                       score.cutoff = 1,
                       KL = "min.exp",
                       quantile.cutoff = 4,
                       BlockOverlap = 0.50,
                       Extension = 0.90,
                       url.link = "https://figshare.com/ndownloader/files/38644505" ,
                       submod.step = 30,
                       min.eRegs = 100
                       ) {

  # Check parameters
  if (!dir.exists(out.dir)) {
    message ("Creating the directory: ", out.dir, " to save the intermediate or final results ...")
    dir.create(out.dir)
  }


  # Quality control
  peaks.use <- rownames(obj[[peak.assay]])[as.vector(as.character(BSgenome::seqnames(obj[[peak.assay]]@ranges)) %in%
                                                       GenomeInfoDb::standardChromosomes(obj[[peak.assay]]@ranges))]
  obj <- subset(obj,
                features = c(
                  rownames(obj[["RNA"]])[!grepl("\\.", rownames(obj[["RNA"]]))],
                  peaks.use
                ))
  message ("Filtered out features not a gene symbol or belonging to non-standard chromosomes,\n",
           "leading to ", nrow(obj[["RNA"]]), " genes and ", nrow(obj[[peak.assay]]), " enhancers.")
  qs::qsave(obj, paste0(out.dir, "Obj_filtered.qsave"))
  message ("Saved the Seurat object after filtering to file: ", paste0(out.dir, "Obj_filtered.qsave"))


  # Libraries
  set.seed(1234)
  if (org == "mm10") {
    require(BSgenome.Mmusculus.UCSC.mm10)
    org.gs <- BSgenome.Mmusculus.UCSC.mm10
  } else if (org == "mm9") {
    require(BSgenome.Mmusculus.UCSC.mm9)
    org.gs <- BSgenome.Mmusculus.UCSC.mm9
  } else if (org == "hg19") {
    require(BSgenome.Hsapiens.UCSC.hg19)
    org.gs <- BSgenome.Hsapiens.UCSC.hg19
  } else {
    require(BSgenome.Hsapiens.UCSC.hg38)
    org.gs <- BSgenome.Hsapiens.UCSC.hg38
  }
  message ("Loaded full genome sequences for ", org, ".")


  # Filter the genes that have annotations
  Seurat::DefaultAssay(obj) <- peak.assay
  obj <- Signac::RegionStats(object = obj, assay = peak.assay,
                     genome = org.gs)
  message ("Computed  the GC content, region lengths, and dinucleotide base frequencies",
           " for regions in the ", peak.assay, " assay and add them to the feature metadata.")
  qs::qsave(obj, paste0(out.dir, "Obj_quality_ctr.qsave"))
  links.df <- filter_nearby_genes(obj = obj, peak.assay = peak.assay)
  obj <- subset(x = obj, features = c(unique(links.df$peak),
                                      unique(links.df$gene)))
  obj <- obj[, Matrix::colSums(obj[['RNA']]@counts) > 0 &
               Matrix::colSums(obj[[peak.assay]]@counts) > 0]
  message ("Removed the cells with zero count of RNA or ATAC, \n",
           "leading to ", ncol(obj), " cells.")
  message ("Filtered the enhancers that are beyond 500 kb from any genes in the Seurat object.\n",
           length(unique(links.df$peak)), " enhancers and ", length(unique(links.df$gene)),
           " genes were retained.")


  # Find highly variable genes and top-ranked enhancers
  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::SCTransform(obj, verbose = T, variable.features.n = var.genes)
  Seurat::DefaultAssay(obj) <- peak.assay
  obj <- Signac::RunTFIDF(obj) # frequency-inverse document frequency (TF-IDF) normalization
  obj <- Signac::FindTopFeatures(obj, min.cutoff = "q5")
  obj <- subset(x = obj, features = c(obj[[peak.assay]]@var.features,
                                        obj[['SCT']]@var.features))
  message ("Saved the Seurat object composed of ", nrow(obj[['RNA']]), " genes (",
           length(obj[["SCT"]]@var.features), " variable genes) and ",
           nrow(obj[[peak.assay]]), " enhancers (", length(obj[[peak.assay]]@var.features),
           " top-ranked enhancers) to file: ", out.dir,
           "Obj_var_genes_top_enhs.qsave")
  qs::qsave(obj, paste0(out.dir, "Obj_var_genes_top_enhs.qsave"))


  # Annotate enhancers with TF binding sites
  load(url(url.link))
  TF.CRE.pairs <- find_TFBS(Seurat::GetAssayData(obj, assay = peak.assay, slot = "data"),
                            TFBS.list = TFBS.list, org = org)
  bound.TFs <- TF.CRE.pairs$CRE
  binding.CREs <- TF.CRE.pairs$TF
  rm(TF.CRE.pairs)
  message ("Identified ", length(binding.CREs), " TFs binding ",
           length(bound.TFs), " enhancers,\n",
           "which were saved to file: ", out.dir, "TF_binding_sites_on_enhs.qsave.")


  # LTMG modeling
  LTMG.obj <- call_LTMG(obj = obj)
  message ("Finished LTMG modeling for ", nrow(LTMG.obj@LTMG@LTMG_discrete),
           " genes across ", ncol(LTMG.obj@LTMG@LTMG_discrete), " cells.\n",
           "There are ", length(unique(range(LTMG.obj@LTMG@LTMG_discrete))),
           " transcriptional regulatory states (TRSs) identified.")


  # Please uncomment the following command and make sure to set a correct working directory
  # so that the following command will generate intermediate files.
  LTMG.obj <- IRISFGM::CalBinaryMultiSignal(LTMG.obj) # partition each gene into TRSs
  LTMG.obj <- IRISFGM::RunBicluster(LTMG.obj, DiscretizationModel = "LTMG", OpenDual = FALSE,
                         NumBlockOutput = n.blocks, BlockOverlap = BlockOverlap,
                         Extension = Extension,
                         BlockCellMin = min.cells)
  message ("Identified ", max(LTMG.obj@BiCluster@CoReg_gene$Condition),
           " IRIS-FGM biclusters from the LTMG matrix.\n",
           "Saved the LTMG object to file: ", out.dir, "LTMG.qsave.")
  qs::qsave(LTMG.obj, paste0(out.dir, "LTMG.qsave"))


  # Get the list of Seurat objects
  rna.dis <- subset(x = LTMG.obj@LTMG@LTMG_discrete,
                    rownames(LTMG.obj@LTMG@LTMG_discrete) %in%
                      rownames(Seurat::GetAssayData(obj, assay = "RNA", slot = "data")))
  atac.dis <- binarize(Seurat::GetAssayData(object = obj, slot = 'data',
                                    assay = peak.assay))
  link.ratio <- sapply(apply(links.df, 2, unique), length)
  obj.list <- subset_object(LTMG.obj = LTMG.obj, object = obj, links.df = links.df,
                            atac.dis = atac.dis, n.blocks = n.blocks,
                            max.peaks = round(nrow(rna.dis) * link.ratio[1] / link.ratio[2]),
                            peak.assay = peak.assay)
  block.list <- split(LTMG.obj@BiCluster@CoReg_gene,
                      f = LTMG.obj@BiCluster@CoReg_gene$Condition) %>%
    sapply(., "[[", "Gene")
  flags <- sapply(obj.list, is.not.null)
  obj.list <- obj.list[flags]
  block.list <- block.list[flags]


  # Construct heterogeneous graphs
  G.list <- build_graph(obj.list = obj.list, obj = obj, rna.dis = rna.dis,
                        atac.dis = atac.dis, ifWeighted = ifWeighted,
                      distance = distance, cicero.covar = cicero.covar,
                      org.gs = org.gs, peak.assay = peak.assay,
                      signac.score = signac.score, signac.pval = signac.pval,
                      min.cells = min.cells)
  flags <- sapply(G.list, is.not.null) # whether the graph is empty
  G.list <- G.list[flags] # filter the graphs
  obj.list <- obj.list[flags] # filter the Seurat objects
  block.list <- block.list[flags]
  qs::qsave(obj.list, paste0(out.dir, "Obj_list.qsave"))
  qs::qsave(G.list, paste0(out.dir, "Graph_list.qsave"))
  if (length(G.list) < 1) {
    stop ("No heterogeneous graph was constructed!")
  }
  message ("List containing ", length(obj.list), " Seurat objects were saved to file: ",
           out.dir, "Obj_list.qsave.")
  message ("Constructed ", length(G.list), " heterogeneous graphs,\n",
           "which were saved to file: ", out.dir, "Graph_list.qsave.")


  # Discover cell-subpopulation-active TF-target pairs
  TFGene.pairs <- find_TF_gene(G.list = G.list, bound.TFs = bound.TFs,
                               binding.CREs = binding.CREs)
  qs::qsave(TFGene.pairs, paste0(out.dir, "TF_target_pairs.qsave"))
  message ("Saved the TF-target pairs to file: ", out.dir, "TF_target_pairs.qsave.")


  # Seeding based upon Steiner forest problem (SFP) model
  seeds <- SFP_seeding(G.list = G.list, obj.list = obj.list, block.list = block.list,
                       bound.TFs = bound.TFs, binding.CREs = binding.CREs,
                       TFGene.pairs = TFGene.pairs, score.cutoff = score.cutoff,
                       quantile.cutoff = quantile.cutoff,
                       rna.dis = rna.dis, atac.dis = atac.dis, KL = KL)
  if (length(seeds) < 1) {
    stop ("No seeds was identified!")
  }
  message (length(seeds), " seeds are identified for hybrid biclustering,\n",
           "which were saved to file: ", out.dir, "Seeds.qsave")
  qs::qsave(seeds, paste0(out.dir, "Seeds.qsave"))


  # Get the list of RNA and ATAC matrices
  rna.list <- get_matrix_list(m = rna.dis, obj.list = obj.list, assay = "RNA") # RNA matrix
  atac.list <- get_matrix_list(m = atac.dis, obj.list = obj.list, assay = peak.assay) # ATAC matrix
  rm(obj.list)


  # Hybrid biclustering
  HBCs <- hybrid_biclust(seeds = seeds, rna.list = rna.list, atac.list = atac.list,
                         top.ngenes = top.ngenes, bound.TFs = bound.TFs,
                         binding.CREs = binding.CREs, G.list = G.list,
                         TFGene.pairs = TFGene.pairs, peak.cutoff = peak.cutoff,
                         c.cutoff = c.cutoff, KL = KL, org.gs = org.gs,
                         intra.cutoff = intra.cutoff, inter.cutoff = inter.cutoff,
                         quantile.cutoff = quantile.cutoff,
                         rna.dis = rna.dis, atac.dis = atac.dis, min.cells = min.cells)
  gene.range <- sapply(HBCs, "[[", "genes") %>% sapply(., length) %>% range
  peak.range <- sapply(HBCs, "[[", "peaks") %>% sapply(., length) %>% range
  cell.range <- sapply(HBCs, "[[", "cells") %>% sapply(., length) %>% range
  message ("Identified ", length(HBCs), " hybrid biclusters (HBCs),\n",
           "which were saved to file: ", out.dir, "HBCs.qsave.\n",
           "HBCs contain ", gene.range[1], "-", gene.range[2], " genes, ",
           peak.range[1], "-", peak.range[2], " enhancers, ",
           cell.range[1], "-", cell.range[2], " cells, and ",
           length(unique(sapply(HBCs, "[[", "TF"))), " TFs.")
  qs::qsave(HBCs, paste0(out.dir, "HBCs.qsave"))


  # Merge significantly overlapped HBCs
  # merged.HBCs <- merge_HBCs(HBCs = HBCs, rna.dis = rna.dis, atac.dis = atac.dis)
  # message (length(merged.HBCs), " HBCs are discovered after fine-tuning.\n")
  # qs::qsave(merged.HBCs, paste0(out.dir, "HBCs_refined.qsave"))


  # Free some memory
  rm(TFGene.pairs)
  rm(atac.list)
  rm(atac.dis)
  rm(HBCs)
  rm(bound.TFs)


  # Optimize the HBCs before submodular optimization
  # obj <- qs::qread(paste0(out.dir, "Obj_var_genes_top_enhs.qsave"))
  # patched.HBCs <- patch_HBCs(merged.HBCs = merged.HBCs, binding.CREs = binding.CREs,
  #                            x = obj, peak.assay = peak.assay,
  #                            distance = patch.dist)
  # qs::qsave(patched.HBCs, paste0(out.dir, "HBCs_optimized.qsave"))


  # Submodular optimization
  if (length(HBCs) <= min.eRegs) {
    submod.HBCs <- HBCs
  } else {
    message("Performing submodular optimization ...")
    sim.m <- compute_sim(HBCs = HBCs) # calculate the pairwise similarity between HBCs
    submod.obj <- sub_mod(HBCs = HBCs, sim.m = sim.m, rna.list = rna.list,
                          G.list = G.list, links.df = links.df,
                          block.list = block.list, n.cells = ncol(rna.dis),
                          obj = obj, min.eRegs = min.eRegs, step = submod.step,
                          peak.assay = peak.assay) # submodular optimization
    rm(sim.m)
    submod.HBCs <- submod.obj$eRegs
    qs::qsave(submod.obj$obj, paste0(out.dir, "Submodular_scores.qsave"))


    return(submod.obj$eRegs)
  }


  # # Extension of HBCs
  # expanded.eGRNs <- expand_eGRNs(obj = obj, submod.HBCs = submod.HBCs,
  #                                peak.assay = peak.assay,
  #                                distance = expand.dist,
  #                                expand.cutoff = expand.cutoff)
  # message ("Discovered ", length(expanded.eGRNs), " enhancer regulons (eRegulons),\n",
  #          "which were saved to file: ", out.dir, "enhancer_regulons.qsave.\n")
  # qs::qsave(expanded.eGRNs, paste0(out.dir, "enhancer_regulons.qsave"))
  #
  #
  # # Return
  # expanded.eGRNs


  submod.HBCs
}



#' Get a simulated jointly profiled scRNA-seq and scATAC-seq data in which
#' several enhancer regulons (eRegulons)
#' are contained
#'
#' @importFrom dplyr %>%
#'
#' @export
#' @rdname create_rna_atac
#'
#' @param obj a \code{Seurat} object used as the prototype to generate simulated dataset
#' other is \code{GRanges} list
#' @param ntfs the number of eRegulons (TFs) to include in the simulated dataset
#' @param ngenes the average number of genes in an eRegulon
#' @param ncells the average number of cells in an eRegulon
#' @param org the organism, e.g., hg38
#' @param atac.assay the assay to save the scATAC-seq data
#' @param gene.links the average number of enhancers linked to each gene in an eRegulon
#' @param distance the maximum distance between a gene and its linked enhancers
#' @param all.genes the number of genes in the simulated \code{Seurat} object
#' @param all.enhs the number of enhancers in the simulated \code{Seurat} object
#' @param all.cells the number of cells in the simulated \code{Seurat} object
#'
#' @return returns a list composed of a eRegulon list and a \code{Seurat} object
create_rna_atac <- function(obj = NULL, ntfs = 5, ngenes = 100,
                            ncells = 100, all.genes = 1000, all.enhs = 3000, all.cells = 1000,
                            org = "hg38", atac.assay = "ATAC", gene.links = 2,
                            distance = 250000, url.link = "https://figshare.com/ndownloader/files/38644505"
                            ) {

  # Parameters
  load(url(url.link))
  ifelse (grepl("^mm", org), sites <- TFBS.list$Mouse, sites <- TFBS.list$Human)
  message ("There are ", length(unique(sites$TF)), " TFs from the JASPAR database.\n",
           ntfs, " eRegulons regulated by different TFs will be generated.\n",
           "The Seurat object contains ", nrow(obj[["RNA"]]),
           " genes, ",
           nrow(obj[[atac.assay]]), " enhancers, and ", ncol(obj), " cells.")


  # Select TFs
  set.seed(123)
  tf.lst <- unique(sites$TF)[sample(seq_along(unique(sites$TF)), size = ntfs)]


  # Select target genes
  # require(dorothea)
  if (grepl("^mm", org)) {
    tf.target <- dorothea::dorothea_mm[, c("tf", "target")]
  } else {
    tf.target <- dorothea::dorothea_hs[, c("tf", "target")]
  }
  message ("Loaded ", length(unique(tf.target$tf)), " TFs collected in DoRothEA database.")
  # require(dplyr)
  tf.target.overlap <- tf.target[tf.target$tf %in% tf.lst,] %>% split(., .$tf)
  message ("There are ", length(names(tf.target.overlap)),
           " common TFs between JASPAR and DoRothEA.")
  ngene.lst <- stats::rpois(ntfs, ngenes)
  gene.lst <- pbmcapply::pbmclapply(seq_along(tf.target.overlap), function(i) {
    x <- tf.target.overlap[[i]]
    x$target[sample(1:nrow(x), min(ngene.lst[i], nrow(x)))]
  }, mc.cores = parallel::detectCores())
  names(gene.lst) <- names(tf.target.overlap)
  message ("Generated regulons composed of ", range(sapply(gene.lst, length))[1],
           "-", range(sapply(gene.lst, length))[2], " genes.")


  # Select enhancers bound by the selected TFs
  # source(paste0(source.dir, "cistrome_tools.R"))
  cis.links <- link_peaks_to_genes(peak.obj = rownames(obj[[atac.assay]]),
                                   gene.obj = Reduce("union", gene.lst), distance = distance,
                                   org = org)
  # require(Signac)
  message ("Identified ", length(unique(Signac::GRangesToString(cis.links))), " enhancers within ",
           distance, " from the selected genes.")
  enhs <- Signac::GRangesToString(cis.links)
  tf.sites <- data.frame(tf = sites$TF[sites$TF %in% tf.lst],
                         peak = Signac::GRangesToString(sites$peak[sites$TF %in% tf.lst])) %>%
    split(., .$tf) %>% sapply(., "[[", 2)


  # Convert TF binding sites to enhancers
  # require(Matrix)
  tf.enhs <- pbmcapply::pbmclapply(tf.sites, function(x) {
    hit.m <- overlap_peak_lst(lst1 = Signac::StringToGRanges(x), lst2 = cis.links)
    # rownames(hit.m) <- StringToGRanges(x)
    # colnames(hit.m) <- enhs
    summ <- Matrix::summary(hit.m)
    enhs[unique(summ$j)]
  }, mc.cores = parallel::detectCores())


  # Ensure each gene was linked to at least one enhancer
  en.regs <- pbmcapply::pbmclapply(names(tf.enhs), function(x) {
    Reduce("c", lapply(gene.lst[[x]], function(y) {
      cis.links[cis.links$gene == y & enhs %in% tf.enhs[[x]]] %>%
        head(n = max(1, rpois(1, gene.links)))
    }))
  }, mc.cores = parallel::detectCores())
  names(en.regs) <- names(tf.enhs)
  message ("Generated ", length(en.regs), " eRegulons composed of ",
           range(sapply(en.regs, length))[1], "-", range(sapply(en.regs, length))[2],
           " enhancer-gene relations.")


  # Build hybrid bicluster (HBCs) based on eRegulons
  rna.m <- Seurat::GetAssayData(object = obj, slot = "counts", assay = "RNA")
  atac.m <- Seurat::GetAssayData(object = obj, slot = "counts", assay = atac.assay)
  hbc.lst <- c()
  ncell.lst <- stats::rpois(length(en.regs), ncells)


  # Initializes the progress bar
  # require(parallel)
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(en.regs), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar


  for (i in seq_along(en.regs)) {
    Sys.sleep(0.1) # Remove this line and add your code


    x <- list(
      TF = names(en.regs)[i],
      genes = unique(en.regs[[i]]$gene),
      enhancers = unique(Signac::GRangesToString(en.regs[[i]]))
    )


    xcells <- sample(1:ncol(obj), ncell.lst[[i]], replace = F)
    rna.m[x$genes, xcells] <- 1 + Reduce("rbind", parallel::mclapply(apply(rna.m[x$genes, xcells], 1, max),
                                                           function(x) rep(x, length(xcells)),
                                                           mc.cores = parallel::detectCores()))
    atac.m[x$enhancers, xcells] <- 1 + Reduce("rbind", parallel::mclapply(apply(atac.m[x$enhancers,
                                                                                       xcells], 1, max),
                                                                function(x) rep(x, length(xcells)),
                                                                mc.cores = parallel::detectCores()))
    x$cells <- colnames(obj)[xcells]
    hbc.lst[[i]] <- x


    utils::setTxtProgressBar(pb, i)
  }


  # Add genes, enhancers, and cells not associated with eRegulons
  hbc.rna <- rna.m[Reduce("union", sapply(hbc.lst, "[[", "genes")),
                   Reduce("union", sapply(hbc.lst, "[[", "cells"))]
  hbc.atac <- atac.m[Reduce("union", sapply(hbc.lst, "[[", "enhancers")),
                   Reduce("union", sapply(hbc.lst, "[[", "cells"))]
  message ("Simulated hybrid biclusters (HBCs) cover ", nrow(hbc.rna),
           " genes, ", nrow(hbc.atac), " enhancers, and ", ncol(hbc.rna),
           " cells.")


  # Add background genes, enhancers, and cells
  choice.genes <- setdiff(rownames(rna.m), rownames(hbc.rna))
  choice.genes <- choice.genes[!grepl("^A", choice.genes) & !grepl("\\.", choice.genes)]
  bg.genes <- sample(choice.genes,
                     all.genes - nrow(hbc.rna))
  bg.enhs <- sample(setdiff(rownames(atac.m), rownames(hbc.atac)),
                    all.enhs - nrow(hbc.atac))
  bg.cells <- sample(setdiff(colnames(atac.m), colnames(hbc.rna)),
                     all.cells - ncol(hbc.rna))
  message ("We have ", length(bg.genes), " genes, ", length(bg.enhs), " enhancers, and ",
           length(bg.cells), " cells as background.")


  # Generate the simulated matrices
  simul.rna <- rna.m[c(rownames(hbc.rna), bg.genes), c(colnames(hbc.rna), bg.cells)]
  simul.atac <- atac.m[c(rownames(hbc.atac), bg.enhs), c(colnames(hbc.atac), bg.cells)]
  rand.rna <- simul.rna[sample(nrow(simul.rna)), sample(ncol(simul.rna))]
  rand.atac <- simul.atac[sample(nrow(simul.atac)), sample(ncol(simul.atac))]
  message ("Generated simulated scRNA-seq and scATAC-seq matrices of sizes: ",
           nrow(rand.rna), " x ", ncol(rand.rna), " and ",
           nrow(rand.atac), " x ", ncol(rand.atac), ".")


  # # Obtain Ensembl based annotation based on organism
  # ensdb <- org_to_DB(org = org)


  # Create the Seurat object
  simul.obj <- rna_atac_matrices_to_Seurat(rna_counts = rand.rna,
                                           atac_counts = rand.atac,
                                           org = org)


  message ("The generated simulated scRNA-seq and scATAC-seq data contains:\n",
           "1. Hybrid biclusters (HBCs) in nested list,\n",
           "2. Seurat object containing the scRNA-seq and scATAC-seq data.")
  list(HBCs = hbc.lst, Seurat = simul.obj)
}



#' Predict cell-type-specific enhancer-driven gene regulatory networks (eGRNs)
#'
#' @importFrom dplyr %>%
#' @importFrom Seurat DefaultAssay SCTransform RunPCA RunUMAP FindMultiModalNeighbors FindClusters
#' @importFrom Signac RunTFIDF FindTopFeatures RunSVD
#'
#' @export
#' @rdname get_cts_en_GRNs
#'
#' @param obj An \code{Seurat} object used to get the eRegulons, NULL by default
#' @param celltype The metadata column indicating the cell types or clusters, NULL by default
#' @param en.regs A list of eRegulons, NULL by default
#' @param peak.assay The chromatin accessibility assay, "ATAC" by default
#' @param rna.dims The number of dimensions for RNA dimension reduction, 50 by default
#' @param atac.dims The number of dimensions for ATAC dimension reduction, 50 by default
#' @param out.dir The directory to save the intermediate results or final results, "./" by default
#' @param padj.cutoff The cutoff of adjusted p-value of hyper-geometric test, 0.05 by default
#' @return Returns a list of eGRNs saved in \code{GRanges} object
#'
get_cts_en_GRNs <- function(obj = NULL, celltype = NULL,
                            en.regs = NULL, peak.assay = "ATAC",
                            rna.dims = 50, atac.dims = 50,
                            padj.cutoff = 0.05,
                            out.dir = "./") {

  # Information
  message ("The Seurat object contains ", nrow(obj[['RNA']]), " genes, ", nrow(obj[[peak.assay]]),
           " enhancers, and ", ncol(obj), " cells.")


  if (!celltype %in% colnames(obj@meta.data)) {
    # RNA analysis
    if (!"SCT" %in% names(obj@assays)) {
      DefaultAssay(obj) <- "RNA"
      obj <- SCTransform(obj, verbose = FALSE)
    } else {
      DefaultAssay(obj) <- "SCT"
    }
    obj <- obj %>% RunPCA() %>% RunUMAP(dims = 1:rna.dims, reduction.name = 'umap.rna',
                                        reduction.key = 'rnaUMAP_')
    message ("Reduced dimensions for transcriptome.")


    # ATAC analysis
    DefaultAssay(obj) <- peak.assay
    if (is.null(obj[[peak.assay]]@var.features)) {
      obj <- RunTFIDF(obj)
      obj <- FindTopFeatures(obj, min.cutoff = 'q0')
    }
    obj <- RunSVD(obj)
    obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:atac.dims,
                   reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    message ("Reduced dimensions for chromatin accessibility.")


    # WNN analysis
    obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                                   dims.list = list(1:rna.dims, 2:atac.dims))
    obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                   reduction.key = "wnnUMAP_")
    obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
    message ("Identified ", length(levels(obj$seurat_clusters)), " cell clusters.\n",
             "Cell clusters are saved in obj$seurat_clusters.\n",
             "The Seurat object after cell clustering is saved in: ", out.dir,
             "Obj_clustered.qsave.")
    celltype <- "seurat_clusters"
  }


  # Hyper-geometric test
  type.cells <- split(obj[[celltype]], f = obj[[celltype]]) %>% sapply(., rownames)
  pval.m <- Reduce("rbind", pbmcapply::pbmclapply(seq_along(en.regs), function(i) {
    group1 <- en.regs[[i]]$cells
    sapply(seq_along(type.cells), function(j) {
      return(phyper(length(intersect(group1, type.cells[[j]])),
             length(type.cells[[j]]),
             ncol(obj) - length(type.cells[[j]]),
             length(group1),
             lower.tail = TRUE))
    })
  }, mc.cores = parallel::detectCores()))
  rownames(pval.m) <- seq_along(en.regs)
  colnames(pval.m) <- names(type.cells)
  message ("Finished hyper-geometric test of ", length(en.regs),
           " eRegulons against ", length(type.cells),
           ' cell types, obtaining p-values ranging from ',
           range(pval.m)[1], " to ", range(pval.m)[2])


  # Bonferroni correction
  padj.m <- pval.m * length(en.regs) * length(type.cells)
  padj.m[padj.m > 1] <- 1.0
  message ("Finished Bonferroni correction and obtained adjusted p-values ranging from ",
           range(padj.m)[1], " to ", range(padj.m)[2])


  # Identify cell-type-specific eRegulons
  cts.reg.ids <- pbmcapply::pbmclapply(1:ncol(padj.m), function(j) {
    return(which(padj.m[, j] < padj.cutoff))
  }, mc.cores = parallel::detectCores())
  names(cts.reg.ids) <- colnames(padj.m)


  # Construct cell-type-specific eRegulons
  cts.reg.ids <- cts.reg.ids[sapply(cts.reg.ids, length) > 0]
  cts.en.grns <- pbmcapply::pbmclapply(names(cts.reg.ids), function(i) {
    links <- unlist(as(lapply(en.regs[cts.reg.ids[[i]]], function(y){
      GenomicRanges::mcols(y$links)$TF <- y$TF
      y$links
    }), "GRangesList"))
    return(list(links = links,
                cell.type = i,
                cells = type.cells[[i]]))
  }, mc.cores = parallel::detectCores())


  return(cts.en.grns)
}



#' Identify cell-type-specific enhancer regulons (eRegulons)
#'
#' @importFrom dplyr %>%
#'
#' @export
#' @rdname get_cts_en_regs
#'
#' @param obj An \code{Seurat} object used to get the eRegulons, NULL by default
#' @param celltype The metadata column indicating the cell types or clusters, NULL by default
#' @param cts.en.grns Cell-type-specific eGRNs, NULL by default
#' @param peak.assay The chromatin accessibility assay, "ATAC" by default
#' @param de.genes A list of differentially expressed genes (DEGs), NULL by default
#' @param out.dir The directory to save the intermediate results or final results, "./" by default
#' @param accessibility whether perform differential accessibility analysis, FALSE by default
#' @param padj.cutoff The cutoff of adjusted p-values of differential expression, 0.05 by default
#'
#' @return Returns a list of cell-type-specific eRegulons and plot dotplot-heatmap
#'
get_cts_en_regs <- function(obj = NULL, peak.assay = "ATAC", de.genes = NULL,
                            cts.en.grns = NULL, accessibility = FALSE, out.dir = "./",
                            min.pct = 0.25, logfc.threshold = 0.25, padj.cutoff = 0.05) {

  # Information
  message ("The Seurat object contains ", nrow(obj[['RNA']]), " genes, ", nrow(obj[[peak.assay]]),
           " enhancers, and ", ncol(obj), " cells.")


  # Differential analyses
  if (is.null(de.genes)) {
    Seurat::DefaultAssay(obj) <- "RNA"
    Idents(obj) <- stats::setNames(obj@meta.data[, celltype], colnames(obj))
    de.genes <- Seurat::FindAllMarkers(obj, only.pos = TRUE,
                                       min.pct = min.pct, logfc.threshold = logfc.threshold)
  }
  de.genes <- de.genes[de.genes$p_val_adj < padj.cutoff,]
  message ("We have ", nrow(de.genes), " differentially expressed genes (DEGs).")
  if (accessibility) {
    Seurat::DefaultAssay(obj) <- peak.assay
    Idents(obj) <- stats::setNames(obj@meta.data[, celltype], colnames(obj))
    da.peaks <- Seurat::FindAllMarkers(obj, only.pos = TRUE,
                                    min.pct = min.pct, logfc.threshold = logfc.threshold)
    da.peaks <- da.peaks[da.peaks$p_val_adj < padj.cutoff,]
    message ("We have ", nrow(da.peaks), " differentially accessible regions (DARs).")
  }


  # Build cell-type-specific eRegulons
  cts.en.regs <- do.call("c", pbmcapply::pbmclapply(cts.en.grns, function(x) {
    links <- x$links
    splitted <- split(links, f = links$TF)
    lapply(seq_along(splitted), function(i) {
      cts.en.reg <- list(
        TF = names(splitted[i]),
        celltype = x$cell.type,
        cells = x$cells,
        genes = intersect(unique(splitted[[i]]$gene),
                          de.genes[de.genes$cluster == x$cell.type, "gene"])
        )
      if (accessibility) {
        cts.en.reg$enhancers <- intersect(unique(Signac::GRangesToString(splitted[[i]])),
                                          da.peaks[da.peaks$cluster == x$cell.type, "gene"])
      } else {
        cts.en.reg$enhancers <- unique(Signac::GRangesToString(splitted[[i]]))
      }
      cts.en.reg$links <- splitted[[i]][Signac::GRangesToString(splitted[[i]]) %in%
                                          cts.en.reg$enhancers &
                                          splitted[[i]]$gene %in% cts.en.reg$genes]
      cts.en.reg$enhancers <- Signac::GRangesToString(cts.en.reg$links)


      return(cts.en.reg)
    })
  }, mc.cores = parallel::detectCores()))
  message ("Identified ", length(cts.en.regs), " cell-type-specific eRegulons.")


  # Visualize cell-type-specific eRegulons using dotplot-heatmap plots
}

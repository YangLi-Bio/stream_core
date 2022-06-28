# To-do list:
# 1. Change the intermediate outputs into text or csv files


# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"


# Source the R scripts
source(paste0(code.dir, "cistrome_tools.R"))
source(paste0(code.dir, "reguome_tools.R"))
source(paste0(code.dir, "epigenome_tools.R"))
source(paste0(code.dir, "multiome_tools.R"))
source(paste0(code.dir, "transcriptome_tools.R"))
source(paste0(code.dir, "utils.R"))
source(paste0(code.dir, "TFBS_list.R"))


# Macros
TOP_TFS <- Inf


# Run STREAM
#' @import Seurat EnsDb.Hsapiens.v86 EnsDb.Mmusculus.v75 BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg38 BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm9 org.Hs.eg.db org.Mm.eg.db Matrix IRISFGM
#' @export
run_stream <- function(obj, var.genes = 3000, top.peaks = 3000,
                       min.cells = 10, out.dir = NULL, org = "hg38",
                       top.ngenes = 15, c.cutoff = 1.0, n.blocks = 500,
                       seed.ratio = 0.30, cicero.covar = 0.00,
                       signac.score = 0.00, min.eGRNs = 100,
                       peak.assay = "ATAC", sim.mode = "both",
                       distance = 500000, signac.pval = 1.0,
                       intra.cutoff = 0.75, inter.cutoff = 0.50,
                       peak.cutoff = 0.5, patch.dist = Inf,
                       cover.blocks = 10, KL = 6,
                       expand.dist = 250000,
                       expand.cutoff = 0.70) {


  # Check parameters
  if (!dir.exists(out.dir)) {
    message ("Creating the directory: ", out.dir, " ..\n")
    dir.create(out.dir)
  }


  # Libraries
  set.seed(1234)
  library(IRISFGM)
  library(Seurat)
  library(Signac)
  if (org == "mm10") {
    library(BSgenome.Mmusculus.UCSC.mm10)
    org.gs <- BSgenome.Mmusculus.UCSC.mm10
  } else if (org == "mm9") {
    library(BSgenome.Mmusculus.UCSC.mm9)
    org.gs <- BSgenome.Mmusculus.UCSC.mm9
  } else if (org == "hg19") {
    library(BSgenome.Hsapiens.UCSC.hg19)
    org.gs <- BSgenome.Hsapiens.UCSC.hg19
  } else {
    library(BSgenome.Hsapiens.UCSC.hg38)
    org.gs <- BSgenome.Hsapiens.UCSC.hg38
  }
  message ("Loading genome sequences of ", org, " ...\n")


  # Filter the genes that have annotations
  DefaultAssay(obj) <- peak.assay
  message ("Computing the GC content, CRE lengths, and dinucleotide base frequencies ...\n")
  obj <- RegionStats(object = obj, assay = peak.assay,
                     genome = org.gs)
  qs::qsave(obj, paste0(out.dir, "Obj_total.qsave"))
  links.df <- filter_nearby_genes(obj = obj, peak.assay = peak.assay)
  # message ("Retained ", nrow(links.df), " CRE-gene pairs which are close to each other.\n")

  obj <- subset(x = obj, features = c(unique(links.df$peak),
                                      unique(links.df$gene)))
  qs::qsave(obj, paste0(out.dir, "Obj_after_filter_nearby_genes.qsave"))
  qs::qsave(links.df, paste0(out.dir, "Nearby_CREs_genes.qsave"))


  # Find highly variable genes and top-ranked CREs
  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(obj, verbose = T, variable.features.n = var.genes)
  DefaultAssay(obj) <- peak.assay
  obj <- RunTFIDF(obj) # frequency-inverse document frequency (TF-IDF) normalization
  obj <- FindTopFeatures(obj, min.cutoff = "q5")
  qs::qsave(obj, paste0(out.dir, "Obj_total.qsave"))
  obj <- subset(x = obj, features = c(obj[[peak.assay]]@var.features,
                                        obj[['SCT']]@var.features))
  message ("Retained ", nrow(obj[['RNA']]), " genes and ", nrow(obj[[peak.assay]]),
           " CREs that are close to at least one CREs/genes.\n")
  qs::qsave(obj, paste0(out.dir, "Obj_varGenes_topCREs.qsave"))


  # Annotate CREs with TF binding sites
  TF.CRE.pairs <- find_TFBS(GetAssayData(obj, assay = peak.assay, slot = "data"),
                            TFBS.list = TFBS.list, org = org)
  qs::qsave(TF.CRE.pairs, paste0(out.dir, "TF_CRE_pairs.qsave"))
  bound.TFs <- TF.CRE.pairs$CRE
  binding.CREs <- TF.CRE.pairs$TF
  rm(TF.CRE.pairs)
  message ("Identified ", length(binding.CREs), " TFs binding ",
           length(bound.TFs), " TFs.\n")


  # LTMG modeling
  LTMG.obj <- call_LTMG(obj = obj)
  LTMG.matrix <- GetLTMGmatrix(LTMG.obj)
  LTMG.matrix <- LTMG.matrix - min(LTMG.matrix)
  qs::qsave(LTMG.matrix, paste0(out.dir, "LTMG_matrix.qsave"))
  qs::qsave(LTMG.obj, paste0(out.dir, "LTMG_obj.qsave"))
  message ("Finished LTMG modeling for ", nrow(LTMG.obj@LTMG@LTMG_discrete),
           " genes across ", ncol(LTMG.obj@LTMG@LTMG_discrete), " cells.\n")
  # p.LTMG <- file.path(tempdir(), "LTMG_matrix.txt")
  # write_LTMG(as.data.frame(binarize(as(LTMG.matrix, "sparseMatrix"))), p.LTMG)


  # Biclustering on transcriptome data
  # setwd(out.dir)
  LTMG.obj <- CalBinaryMultiSignal(LTMG.obj)
  # LTMG.obj@LTMG@LTMG_BinaryMultiSignal <- LTMG.obj@LTMG@LTMG_discrete -
  #   min(LTMG.obj@LTMG@LTMG_discrete)
  block.original <- RunBicluster(
    LTMG.obj,
    DiscretizationModel = "LTMG", OpenDual = F,
    NumBlockOutput = n.blocks, BlockOverlap = 0.25,
    BlockCellMin = min.cells, Extension = c.cutoff)
  rm(LTMG.obj)
  qs::qsave(block.original, paste0(out.dir, "QUBIC_block_obj.qsave"))
  block.list <- load_blocks(block.original = block.original,
                            cover.blocks = cover.blocks, n.blocks = n.blocks,
                            sim.mode = sim.mode, rank.blocks = T)
  qs::qsave(block.list, paste0(out.dir, "QUBIC_blocks.qsave"))
  message ("Identified ", length(block.list),
           " QUBIC biclusters.\n")


  # Get the list of Seurat objects
  rna.dis <- subset(x = LTMG.matrix, rownames(LTMG.matrix) %in%
                      rownames(GetAssayData(obj, assay = "RNA",
                                            slot = "data")))
  atac.dis <- binarize(GetAssayData(object = obj, slot = 'data',
                                    assay = peak.assay))
  obj.list <- subset_object(block.list = block.list, object = obj, links.df = links.df,
                            atac.dis = atac.dis, max.peaks = top.peaks,
                            min.cells = 0, peak.assay = peak.assay)
  # flags <- sapply(obj.list, function(x) {
  #   if (nrow(x[["RNA"]]) < 1 |
  #       nrow(x[[peak.assay]]) < 1) {
  #     return(F)
  #   }
  #   return(T)
  # })
  flags <- sapply(obj.list, is.not.null)
  obj.list <- obj.list[flags]
  block.list <- block.list[flags]
  if (length(block.list) < 1) {
    stop ("No QUBIC biclusters can be used to construct Steiner forest!\n")
  }
  message (length(block.list), " QUBIC biclusters will be used for Steiner forest modeling.\n")


  # Construct heterogeneous graphs
  G.list <- build_graph(obj.list = obj.list,
                      distance = distance, cicero.covar = cicero.covar,
                      org.gs = org.gs,
                      signac.score = signac.score, signac.pval = signac.pval,
                      min.cells = 0)
  flags <- sapply(G.list, is.not.null) # whether the graph is empty
  G.list <- G.list[flags] # filter the graphs
  block.list <- block.list[flags] # filter the blocks
  obj.list <- obj.list[flags] # filter the seurat objects
  qs::qsave(block.list, paste0(out.dir, "QUBIC_biclusters.qsave"))
  qs::qsave(rna.dis, paste0(out.dir, "RNA_discretized.qsave"))
  qs::qsave(atac.dis, paste0(out.dir, "ATAC_binarized.qsave"))
  qs::qsave(obj.list, paste0(out.dir, "Obj_list.qsave"))
  qs::qsave(G.list, paste0(out.dir, "Graph_list.qsave"))
  if (length(G.list) < 1) {
    stop ("No heterogeneous graph is constructed!\n")
  }
  message ("Finished constructing ", length(G.list), " heterogeneous graphs.\n")


  # Discover cell-subpopulation-active TF-target pairs
  TFGene.pairs <- find_TF_gene(G.list = G.list, bound.TFs = bound.TFs,
                               binding.CREs = binding.CREs)
  qs::qsave(TFGene.pairs, paste0(out.dir, "Cell_subpopulation_active_TF_target.qsave"))


  # Seeding based upon Steiner forest problem (SFP) model
  seeds <- SFP_seeding(block.list = block.list, G.list = G.list, obj.list = obj.list,
                       bound.TFs = bound.TFs, binding.CREs = binding.CREs,
                       TFGene.pairs = TFGene.pairs, score.cutoff = 1,
                       rna.dis = rna.dis, atac.dis = atac.dis, KL = KL)
  if (length(seeds) < 1) {
    stop ("No seeds is identified!\n")
  }
  message (length(seeds), " seeds are identified for hybrid biclustering.\n")
  qs::qsave(seeds, paste0(out.dir, "Seeds.qsave"))


  # Get the list of RNA and ATAC matrices
  rna.list <- get_matrix_list(m = rna.dis, obj.list = obj.list, assay = "RNA") # RNA matrix
  atac.list <- get_matrix_list(m = atac.dis, obj.list = obj.list, assay = peak.assay) # ATAC matrix
  rm(obj.list)
  qs::qsave(rna.list, paste0(out.dir, "RNA_matrix_list.qsave"))
  qs::qsave(atac.list, paste0(out.dir, "ATAC_matrix_list.qsave"))


  # Hybrid biclustering
  HBCs <- hybrid_biclust(seeds = seeds, rna.list = rna.list, atac.list = atac.list,
                         top.ngenes = top.ngenes, bound.TFs = bound.TFs,
                         binding.CREs = binding.CREs, G.list = G.list,
                         TFGene.pairs = TFGene.pairs, peak.cutoff = peak.cutoff,
                         c.cutoff = c.cutoff, KL = KL, org.gs = org.gs,
                         intra.cutoff = intra.cutoff, inter.cutoff = inter.cutoff,
                         rna.dis = rna.dis, atac.dis = atac.dis, min.cells = min.cells)
  message ("Identified ", length(HBCs), " hybrid biclusters (HBCs).\n")
  qs::qsave(HBCs, paste0(out.dir, "HBCs.qsave"))


  # Merge significantly overlapped HBCs
  merged.HBCs <- merge_HBCs(HBCs = HBCs, rna.dis = rna.dis, atac.dis = atac.dis)
  message (length(merged.HBCs), " HBCs are discovered after fine-tuning.\n")
  qs::qsave(merged.HBCs, paste0(out.dir, "HBCs_refined.qsave"))


  # Free some memory
  rm(TFGene.pairs)
  rm(atac.list)
  rm(atac.dis)
  rm(HBCs)
  rm(bound.TFs)


  # Optimize the HBCs before submodular optimization
  obj <- qs::qread(paste0(out.dir, "Obj_total.qsave"))
  patched.HBCs <- patch_HBCs(merged.HBCs = merged.HBCs, binding.CREs = binding.CREs,
                             x = obj, peak.assay = peak.assay,
                             distance = patch.dist)
  qs::qsave(patched.HBCs, paste0(out.dir, "HBCs_optimized.qsave"))


  # Submodular optimization
  if (length(patched.HBCs) <= min.eGRNs) {
    submod.HBCs <- patched.HBCs
  } else {
    message("Performing submodular optimization ...\n")
    sim.m <- compute_sim(HBCs = patched.HBCs) # calculate the pairwise similarity between HBCs
    submod.obj <- sub_mod(HBCs = patched.HBCs, sim.m = sim.m, rna.list = rna.list,
                          G.list = G.list,
                          block.list = block.list, n.cells = ncol(rna.dis),
                          obj = obj, min.eGRNs = min.eGRNs,
                          peak.assay = peak.assay, distance = distance) # submodular optimization
    rm(sim.m)
    submod.HBCs <- submod.obj$eGRNs
    qs::qsave(submod.obj$obj, paste0(out.dir, "Submodular_scores.qsave"))
  }
  message ("Submodular optimization identified ", length(submod.HBCs),
           " enhancer gene regulatory networks (eGRNs).\n")
  qs::qsave(submod.HBCs, paste0(out.dir, "HBCs_submod.qsave"))


  # Extension of eGRNs
  expanded.eGRNs <- expand_eGRNs(obj = obj, submod.HBCs = submod.HBCs,
                                 peak.assay = peak.assay,
                                 distance = expand.dist,
                                 expand.cutoff = expand.cutoff)


  # Return
  expanded.eGRNs
}

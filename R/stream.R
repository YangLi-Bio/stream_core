# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"


# Source the R scripts
source(paste0(code.dir, "cistrome_tools.R"))
source(paste0(code.dir, "reguome_tools.R"))
source(paste0(code.dir, "utils.R"))


# Run STREAM
#' @import Seurat EnsDb.Hsapiens.v86 EnsDb.Mmusculus.v75 BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg38 BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm9 org.Hs.eg.db org.Mm.eg.db Matrix IRISFGM
#' @export
run_stream <- function(obj, var.genes = 3000, top.peaks = 3000,
                       min.cells = 10, out.dir = NULL, org = "hg38",
                       top.ngenes = 15, c.cutoff = 1.0, n.blocks = 500,
                       seed.ratio = 0.30, civero.covar = 0.00,
                       signac.score = 0.00, min.eGRNs = 100,
                       peak.assay = "ATAC", sim.mode = "both",
                       cover.blocks = 10, KL = 6) {
  set.seed(1234)


  # Libraries
  library(IRISFGM)
  if (org == "mm10") {
    org.gs <- BSgenome.Mmusculus.UCSC.mm10
    library(BSgenome.Mmusculus.UCSC.mm10)
  } else if (org == "mm9") {
    org.gs <- BSgenome.Mmusculus.UCSC.mm9
    library(BSgenome.Mmusculus.UCSC.mm9)
  } else if (org == "hg19") {
    org.gs <- BSgenome.Hsapiens.UCSC.hg19
    library(BSgenome.Hsapiens.UCSC.hg19)
  } else {
    org.gs <- BSgenome.Hsapiens.UCSC.hg38
    library(BSgenome.Hsapiens.UCSC.hg38)
  }
  message ("Loading genome sequences of ", org, " ...\n")


  # Filter the genes that have annotations
  DefaultAssay(obj) <- peak.assay
  message ("Computing the GC content, CRE lengths, and dinucleotide base frequencies ...\n")
  obj <- RegionStats(object = obj, assay = peak.assay,
                     genome = org.gs)
  links.df <- filter_nearby_genes(obj = obj, peak.assay = peak.assay)
  # message ("Retained ", nrow(links.df), " CRE-gene pairs which are close to each other.\n")

  obj <- subset(x = obj, features = c(unique(links.df$peak),
                                      unique(links.df$gene)))
  qs::qsave(obj, paste0(out.dir, "Obj_after_filter_nearby_genes.qsave"))
  qs::qsave(links.df, paste0(out.dir, "Nerby_CREs_genes.qsave"))


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
                            org = org)
  bound.TFs <- TF.CRE.pairs$CRE
  binding.CREs <- TF.CRE.pairs$TF
  rm(TF.CRE.pairs)
  message ("Identified ", length(binding.CREs), " TFs binding ",
           length(bound.TFs), " TFs.\n")


  # LTMG modeling
  LTMG.matrix <- call_LTMG(obj = obj)
  qs::qsave(LTMG.matrix, paste0(out.dir, "LTMG_matrix.qsave"))
  # p.LTMG <- file.path(tempdir(), "LTMG_matrix.txt")
  # write_LTMG(as.data.frame(binarize(as(LTMG.matrix, "sparseMatrix"))), p.LTMG)


  # Biclustering on transcriptome data
  block.original <- RunBicluster(
    CreateIRISFGMObject(as.data.frame(binarize(as(LTMG.matrix, "sparseMatrix")))),
    DiscretizationModel = "LTMG", OpenDual = F,
    NumBlockOutput = n.blocks, BlockOverlap = 0.25,
    BlockCellMin = min.cells, Extension = c.cutoff)
  message ("Identified ", length(unique(block.list[, 2, drop = F])),
           " QUBIC biclusters.\n")
  block.list <- retain_blocks(block.original = block.original, sim.mode = sim.mode,
                              cover.blocks = cover.blocks)
  # run_QUBIC(path = qubic.path, file = ltmg.file, f = 0.25,
  #   dual = dual.mode, k = min.cells, o = n.blocks)


  # Get the list of Seurat objects
  rna.dis <- subset(x = LTMG.matrix, rownames(LTMG.matrix) %in%
                      rownames(GetAssayData(obj, assay = "RNA",
                                            slot = "data")))
  atac.dis <- binarize(GetAssayData(object = obj, slot = 'data',
                                    assay = peak.assay))
  obj.list <- subset_object(block.list = block.list, object = obj, links.df = links.df,
                            atac.dis = atac.dis, max.peaks = top.peaks,
                            min.cells = 0)
  flags <- sapply(obj.list, is.not.full)
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
  flags <- !unlist(sapply(G.list, is.null)) # whether the graph is empty
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
                       TFGene.pairs = TFGene.pairs,
                       rna.dis = rna.dis, atac.dis = atac.dis, KL = KL)
  if (length(seeds) < 1) {
    stop ("No seeds is identified!\n")
  }
  mesage (length(seeds), " seeds are identified for hybrid biclustering.\n")
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
                         binding.CREs = binding.CREs, G.list = G.list, TFGene.pairs = TFGene.pairs,
                         c.cutoff = c.cutoff, KL = KL, org.gs = org.gs,
                         rna.dis = rna.dis, atac.dis = atac.dis, min.cells = min.cells)
}

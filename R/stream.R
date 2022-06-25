# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"


# Source the r scripts
source(paste0(code.dir, "cistrome_tools.R"))


# Run STREAM
#' @import Seurat EnsDb.Hsapiens.v86 EnsDb.Mmusculus.v75 BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg38 BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm9 org.Hs.eg.db org.Mm.eg.db
#' @export
run_stream <- function(obj, var.genes = 3000, top.peaks = 3000, 
                       min.cells = 10, out.dir = NULL, org = "hg38", 
                       top.ngenes = 15, c.cutoff = 1.0, n.blocks = 500, 
                       seed.ratio = 0.30, civero.covar = 0.00, 
                       signac.score = 0.00, min.eGRNs = 100, 
                       peak.assay = "ATAC") {
  set.seed(1234)
  
  
  # Libraries
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
  
  
  # Filter the genes that have annotations
  DefaultAssay(obj) <- peak.assay
  obj <- RegionStats(object = obj, assay = peak.assay, 
                     genome = org.gs)
  links.df <- filter_nearby_genes(obj = obj, peak.assay = peak.assay)
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
  obj.all <- obj
  obj <- subset(x = obj, features = c(obj[[peak.assay]]@var.features,
                                        obj[['SCT']]@var.features))
  qs::qsave(obj, paste0(out.dir, "Obj_varGenes_topCREs.qsave"))
  
  
  # Annotate CREs with TF binding sites
  
}

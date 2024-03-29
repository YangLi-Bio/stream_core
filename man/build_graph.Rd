% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiome_tools.R
\name{build_graph}
\alias{build_graph}
\title{Build a list of heterogeneous graphs based on the input list of \code{Seurat} objects}
\usage{
build_graph(
  obj.list,
  obj = NULL,
  rna.dis,
  atac.dis,
  distance = 250000,
  cicero.covar = -Inf,
  org.gs = BSgenome.Hsapiens.UCSC.hg38,
  signac.score = -Inf,
  signac.pval = Inf,
  min.cells = 10,
  ifWeighted = F,
  peak.assay = "ATAC"
)
}
\description{
Build a list of heterogeneous graphs based on the input list of \code{Seurat} objects
}
\keyword{internal}

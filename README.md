# STREAM
We present a new algorithm, STREAM, for enhancer-driven gene regulatory network (eGRN) inference from transcriptome and chromatin accessibility profiled from the same single cell populations. The algorithm substantially improves the prediction accuracy of relations among transcription factors (TFs), enhancers, and genes, by achieving global optimization based on two key new ideas: (i) we developed the Steiner forest problem (SFP  ) model based on a heterogeneous graph to identify the set of highly confident enhancer-gene relations which underlie a context-specific functional gene module (FGM); and (ii) we designed a hybrid biclustering pipeline integrated with submodular optimization for inferring eGRNs by identifying the optimal subset from a set of hybrid biclusters (HBCs), each of which represents co-regulated genes by the same TF, and co-accessible enhancers bound by the same TF, occurring over a subset of cells. 

# Workflow of STREAM
STREAM combines the Steiner forest problem (SFP) model and submodular optimization, respectively, to discover the enhancer-gene relations and TF-enhancer-gene relations in a global optimization manner.

![Figure_1_workflow_12122022](https://user-images.githubusercontent.com/35290254/207700357-5bc15019-4733-48b1-ad3f-7613769b8285.jpg)

# Prerequisites
- Seurat
- Signac
- cicero
- qualV
- GenomicRanges
- ensdb.hsapiens.v75
- ensdb.hsapiens.v86
- ensdb.mmusculus.v75
- ensdb.mmusculus.v79
- bsgenome.mmusculus.ucsc.mm9
- bsgenome.mmusculus.ucsc.mm10
- bsgenome.hsapiens.ucsc.hg19
- bsgenome.hsapiens.ucsc.hg38

# Installation
The latest developmental version of STREAM can be downloaded from GitHub and installed from source by `devtools::install_github('YangLi-Bio/stream_core')`.

# License
Scissor is licensed under the GNU General Public License v3.0.

Improvements and new features of STREAM will be updated on a regular basis. Please post on the GitHub discussion page with any questions.

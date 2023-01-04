# link_cor <- function(x.atac, distance = 500000, cicero.covar = 0) {
#
#   sep <- c("-", "-")
#   GR.ll <- Signac::StringToGRanges(rownames(x.atac)) # transform strings to GRanges
#   GR.expand <- GenomicRanges::resize(x = GR.ll, width = distance, fix = 'center') # find the TSS location
#   overlaps <- GenomicAlignments::findOverlaps(
#     query = GR.ll,
#     subject = GR.expand,
#     type = 'any',
#     select = 'all'
#   ) # find the peaks overlaped with the extended genomic ranges of peaks
#   message ("Identified overlapped enhancers within ", distance, " bps.\n")
#   hit_matrix <- Matrix::sparseMatrix(
#     i = queryHits(x = overlaps),
#     j = subjectHits(x = overlaps),
#     x = 1,
#     dims = c(length(x = GR.ll), length(x = GR.expand))
#   ) # build a sparse matrix to record the overlaps between peaks and extended genomic ranges of genes
#   rownames(x = hit_matrix) <- Signac::GRangesToString(grange = GR.ll, sep = sep) # use peak names as the row names
#   colnames(x = hit_matrix) <- Signac::GRangesToString(grange = GR.ll, sep = sep) # use peak names as the column names
#   hit_matrix[!Ryacas::upper.tri(hit_matrix, diag = F)] <- 0 # mask the diagonal and upper triangular elements
#   row.col <- which(hit_matrix == 1, arr.ind = T) # double vector to hold the row and column IDs for elements
#   # equalling one
#   if (nrow(row.col) < 1) {
#     message ('No enhancer-enhancer relations was discovered.')
#     return(NULL)
#   }
#
#
#   # Assign covariance as the matrix elements
#   for (i in 1 : nrow(row.col)) {
#     hit_matrix[row.col[i, 1], row.col[i, 2]] <- cov(x.atac[row.col[i, 1], ],
#                                                     x.atac[row.col[i, 2], ])
#   }
#   hit_matrix[which(hit_matrix < cicero.covar)] <- 0 # discard the negative elements
#   summ <- Matrix::summary(hit_matrix) # convert the matrix into a sparse matrix
#   cor.links <- data.frame(Origin = rownames(hit_matrix)[summ$i],
#                           Destination = colnames(hit_matrix)[summ$j],
#                           Weight      = summ$x) # transform the sparse matrix into a data frame
#   colnames(cor.links) <- c('node1', 'node2', 'weight') # name the columns
#   max.weight <- max(cor.links$weight)
#   min.weight <- min(cor.links$weight)
#   diff <- max.weight - min.weight
#   quiet(ifelse(diff > 0, cor.links$weight <- (max.weight - cor.links$weight) /
#                  diff, cor.links$weight <- 0)) # normalize the weights
#   message ('Identified ', nrow(cor.links), ' enhancer-enhancer relations.')
#
#
#   cor.links
# }



# link_cicero <- function(x, distance = 500000, cicero.covar = 0,
#                         org.gs = BSgenome.Hsapiens.UCSC.hg38) {
#
#   summ <- Matrix::summary(x) # convert the matrix into a sparse matrix
#   cicero.data <- data.frame(Origin = rownames(x)[summ$i],
#                             Destination = colnames(x)[summ$j],
#                             Weight      = summ$x) # transform the sparse matrix into a data frame
#   input.cds <- make_atac_cds(cicero.data, binarize = F) %>% detect_genes
#   input.cds <- input.cds[Matrix::rowSums(exprs(input.cds)) != 0, ] %>% estimate_size_factors %>%
#     preprocess_cds(method = "LSI", verbose = F) %>% reduce_dimension(reduction_method = 'UMAP',
#                                                                      preprocess_method = "LSI")
#   umap.coords <- SingleCellExperiment::reducedDims(input.cds)$UMAP # obtain the UMAP coordinates
#   cicero.cds <- make_cicero_cds(input.cds, reduced_coordinates = umap.coords)
#   genome.info <- data.frame(org.gs@seqinfo@seqnames,
#                             org.gs@seqinfo@seqlengths) # genome sequence lengths
#   colnames(genome.info) <- c("seqnames", "seqlengths") # rename the columns
#   cicero.links <- run_cicero(cicero.cds, genome.info,
#                              window = distance) # build peak-peak linkages using cicero
#   colnames(cicero.links) <- c('node1', 'node2', 'weight')
#   cicero.links$node2 <- as.character(cicero.links$node2) # convert factors into characters
#   cicero.links <- data.table::rbindlist(apply(cicero.links, 1, function(r) {
#     ifelse(r[1] >= r[2], return(list(node1 = r[1], node2 = r[2], weight = r[3])),
#            return(list(node1 = r[2], node2 = r[1], weight = r[3])))
#   }), fill = T) %>% dplyr::distinct()
#
#
#   cicero.links$weight <- as.numeric(cicero.links$weight) # convert characters into numeric values
#   cicero.links$weight[is.na(cicero.links$weight)] <- 0 # substitute the NA values
#   max.weight <- max(cicero.links$weight)
#   min.weight <- min(cicero.links$weight)
#   diff <- max.weight - min.weight
#   ifelse(diff > 0, cicero.links$weight <- (max.weight - cicero.links$weight) /
#            diff, cicero.links$weight <- 0) # normalize the weights
#   message ('Generated ', nrow(cicero.links), 'enhancer-enhancer relations.')
#
#   cicero.links
# }



#' Calculate whether two list of \code{GRanges} are overlapped
#'
#' @keywords internal
#'
overlap_peak_lst <- function(lst1, lst2) {

  overlaps <- GenomicAlignments::findOverlaps(
    query = lst1,
    subject = lst2,
    type = 'any',
    select = 'all'
  ) # find the peaks overlapped with the extended genomic ranges of genes
  # require(Matrix)


  return(Matrix::sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = lst1), length(x = lst2))
  )) # build a sparse matrix to record the overlaps between peaks and extended genomic ranges of genes
}



#' Function org_to_DB used to get the database based on organism
#'
#' @keywords internal
#'
org_to_DB <- function(org = "hg38") {

  message ("The organism version is ", org, ".")


  # Return the database
  if (grepl("mm", org)) {
    message ("Loading the database EnsDb.Mmusculus.v79 ...")
    require(EnsDb.Mmusculus.v79)
    return(EnsDb.Mmusculus.v79)
  } else {
    message ("Loading the database EnsDb.Hsapiens.v86 ...")
    require(EnsDb.Hsapiens.v86)
    return(EnsDb.Hsapiens.v86)
  }
}



#' Generate genomic windows
#'
#' @keywords internal
#'
#' @import dplyr
#'
generate_windows <- function (window, genomic_coords) {

  if (!is(genomic_coords, "data.frame")) {
    chr_maxes <- read.table(genomic_coords)
  }
  else {
    chr_maxes <- genomic_coords
  }
  names(chr_maxes) <- c("V1", "V2")
  win_ranges <- plyr::ddply(chr_maxes, plyr::.(V1), function(x) {
    r <- seq(from = 1, to = x$V2[1], by = window/2)
    l <- r + window - 1
    data.frame(start = r, end = l)
  })
  gr <- GenomicRanges::GRanges(win_ranges$V1,
                               ranges = IRanges::IRanges(win_ranges$start,
                                                         win_ranges$end))


  return(gr)
}



#' Get genomic ranges
#'
#' @keywords internal
#'
get_genomic_range <- function (grs, cds, win) {

  end1 <- as.numeric(as.character(GenomicRanges::end(grs[win])))
  end2 <- as.numeric(as.character(GenomicRanges::start(grs[win])))
  win_range <- cds[(monocle3::fData(cds)$bp1 < end1 & monocle3::fData(cds)$bp1 >
                      end2) | (monocle3::fData(cds)$bp2 < end1 & monocle3::fData(cds)$bp2 > end2),
  ]
  win_range <- win_range[as.character(monocle3::fData(win_range)$chr) ==
                           gsub("chr", "", as.character(GenomicRanges::seqnames(grs[win]))),
  ]
  monocle3::fData(win_range)$mean_bp <- (as.numeric(as.character(monocle3::fData(win_range)$bp1)) +
                                 as.numeric(as.character(monocle3::fData(win_range)$bp2)))/2


  return(win_range)
}



#' Find distance parameter
#'
#' @keywords internal
#'
find_distance_parameter <- function(dist_mat,
                                    gene_range,
                                    maxit,
                                    null_rho,
                                    s,
                                    distance_constraint,
                                    distance_parameter_convergence) {

  if (sum(dist_mat > distance_constraint) / 2 < 1) {
    return("No long edges")
  }

  found <- FALSE
  starting_max <- 2
  distance_parameter <- 2
  distance_parameter_max <- 2
  distance_parameter_min <- 0
  it <- 0
  while(found != TRUE & it < maxit) {
    vals <- monocle3::exprs(gene_range)
    cov_mat <- cov(as.data.frame(t(vals)))
    diag(cov_mat) <- Matrix::diag(cov_mat) + 1e-4

    rho <- get_rho_mat(dist_mat, distance_parameter, s)

    GL <- glasso::glasso(cov_mat, rho)
    big_entries <- sum(dist_mat > distance_constraint)

    if (((sum(GL$wi[dist_mat > distance_constraint] != 0) / big_entries) > 0.05) |
        (sum(GL$wi == 0) / (nrow(GL$wi) ^ 2) < 0.2 ) ) {
      longs_zero <- FALSE
    } else {
      longs_zero <- TRUE
    }

    if (longs_zero != TRUE | (distance_parameter == 0)) {
      distance_parameter_min <- distance_parameter
    } else {
      distance_parameter_max <- distance_parameter
    }
    new_distance_parameter <- (distance_parameter_min +
                                 distance_parameter_max)/2

    if(new_distance_parameter == starting_max) {
      new_distance_parameter <- 2 * starting_max
      starting_max <- new_distance_parameter
    }

    if (distance_parameter_convergence > abs(distance_parameter -
                                             new_distance_parameter)) {
      found <- TRUE
    } else {
      distance_parameter <- new_distance_parameter
    }
    it <- it + 1
  }
  if (maxit == it) warning ("maximum iterations hit")


  return(distance_parameter)
}



#' Get distance matrix
#'
#' @keywords internal
#'
calc_dist_matrix <- function (gene_range) {

  dist_mat <- Matrix::as.matrix(stats::dist(monocle3::fData(gene_range)$mean_bp))
  row.names(dist_mat) <- colnames(dist_mat) <- row.names(monocle3::fData(gene_range))


  return(dist_mat)
}



#' Estimate distance parameter using \code{cicero}
#'
#' @keywords internal
#'
estimate_distance_parameter <- function(cds, window = 5e+05, maxit = 100, s = 0.75, sample_num = 100,
                                         distance_constraint = 250000, distance_parameter_convergence = 1e-22,
                                         max_elements = 200, genomic_coords = cicero::human.hg19.genome,
                                         max_sample_windows = 500) {

  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(assertthat::is.number(window))
  assertthat::assert_that(assertthat::is.count(maxit))
  assertthat::assert_that(assertthat::is.number(s), s < 1,
                          s > 0)
  assertthat::assert_that(assertthat::is.count(sample_num))
  assertthat::assert_that(assertthat::is.count(distance_constraint))
  assertthat::assert_that(distance_constraint < window)
  assertthat::assert_that(assertthat::is.number(distance_parameter_convergence))
  if (!is.data.frame(genomic_coords)) {
    assertthat::is.readable(genomic_coords)
  }
  assertthat::assert_that(assertthat::is.count(max_sample_windows))
  grs <- generate_windows(window, genomic_coords)
  fData(cds)$chr <- gsub("chr", "", fData(cds)$chr)
  fData(cds)$bp1 <- as.numeric(as.character(monocle3::fData(cds)$bp1))
  fData(cds)$bp2 <- as.numeric(as.character(monocle3::fData(cds)$bp2))
  distance_parameters <- list()
  distance_parameters_calced <- 0
  it <- 0
  while (sample_num > distance_parameters_calced &
         it < max_sample_windows) {
    it <- it + 1
    win <- sample(seq_len(length(grs)), 1)
    GL <- "Error"
    win_range <- get_genomic_range(grs, cds, win)
    if (nrow(monocle3::exprs(win_range)) <= 1) {
      (next)()
    }
    if (nrow(monocle3::exprs(win_range)) > max_elements) {
      (next)()
    }
    dist_matrix <- calc_dist_matrix(win_range)
    distance_parameter <- find_distance_parameter(dist_matrix,
                                                  win_range, maxit = maxit,
                                                  null_rho = 0, s,
                                                  distance_constraint = distance_constraint,
                                                  distance_parameter_convergence = distance_parameter_convergence)
    if (!is(distance_parameter, "numeric")) {
      next
    }
    distance_parameters = c(distance_parameters, distance_parameter)
    distance_parameters_calced <- distance_parameters_calced +
      1
  }
  if (length(distance_parameters) < sample_num)
    warning(paste0("Could not calculate sample_num distance_parameters (",
                   length(distance_parameters), " were calculated) - see ",
                   "documentation details"))
  if (length(distance_parameters) == 0)
    stop("No distance_parameters calculated")


  unlist(distance_parameters)
}



#' Get rho mat
#'
#' @keywords internal
#'
get_rho_mat <- function(dist_matrix, distance_parameter, s) {
  xmin <- 1000
  out <- (1-(xmin/dist_matrix)^s) * distance_parameter
  out[!is.finite(out)] <- 0
  out[out < 0] <- 0


  out
}


#' Generate \code{cicero} model
#'
#' @keywords internal
#'
generate_cicero_models <- function (cds, distance_parameter, s = 0.75, window = 5e+05,
                                    max_elements = 200, genomic_coords = NULL)
{
  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(assertthat::is.number(distance_parameter))
  assertthat::assert_that(assertthat::is.number(s), s < 1,
                          s > 0)
  assertthat::assert_that(assertthat::is.number(window))
  assertthat::assert_that(assertthat::is.count(max_elements))
  if (!is.data.frame(genomic_coords)) {
    assertthat::is.readable(genomic_coords)
  }
  grs <- generate_windows(window, genomic_coords)
  fData(cds)$chr <- gsub("chr", "", fData(cds)$chr)
  fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
  fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))
  outlist <- pbmcapply::pbmclapply(seq_len(length(grs)), mc.cores = max(parallel::detectCores() / 2, 1),
                                function(win) {
                                  GL <- "Error"
                                  win_range <- get_genomic_range(grs, cds, win)
                                  if (nrow(exprs(win_range)) <= 1) {
                                    return("Zero or one element in range")
                                  }
                                  if (nrow(exprs(win_range)) > max_elements) {
                                    return("Too many elements in range")
                                  }
                                  dist_matrix <- calc_dist_matrix(win_range)
                                  rho_mat <- get_rho_mat(dist_matrix, distance_parameter,
                                                         s)
                                  vals <- Matrix::as.matrix(exprs(win_range))
                                  cov_mat <- cov(t(vals))
                                  Matrix::diag(cov_mat) <- diag(cov_mat) + 1e-04
                                  GL <- glasso::glasso(cov_mat, rho_mat)
                                  colnames(GL$w) <- row.names(GL$w) <- row.names(vals)
                                  colnames(GL$wi) <- row.names(GL$wi) <- row.names(vals)
                                  return(GL)
                                })
  names_df <- as.data.frame(grs)
  names(outlist) <- paste(names_df$seqnames, names_df$start,
                          names_df$end, sep = "_")


  outlist
}



#' Assemble \code{cicero} coaccessibility relations
#'
#' @keywords internal
#'
assemble_connections <- function (cicero_model_list, silent = FALSE) {

  types <- vapply(cicero_model_list, FUN = class, FUN.VALUE = "character")
  char_hbn <- cicero_model_list[types == "character"]
  gl_only <- cicero_model_list[types == "list"]
  if (!silent) {
    print(paste("Successful cicero models: ", length(gl_only)))
    print("Other models: ")
    print(table(unlist(char_hbn)))
    print(paste("Models with errors: ", sum(is.null(cicero_model_list))))
  }
  cors <- lapply(gl_only, function(gl) {
    cors <- stats::cov2cor(gl$w)
    data.table::melt(as.data.table(cors, keep.rownames = TRUE),
                     measure = patterns("[0-9]"))
  })
  cors <- data.table::rbindlist(cors)
  names(cors) <- c("Var1", "Var2", "value")
  data.table::setkey(cors, "Var1", "Var2")
  cors_rec <- as.data.frame(cors[, list(mean_coaccess = reconcile(value)),
                                 by = "Var1,Var2"])
  names(cors_rec) <- c("Peak1", "Peak2", "coaccess")
  cors_rec <- cors_rec[cors_rec$Peak1 != cors_rec$Peak2, ]


  cors_rec
}



#' Split peak names
#'
#' @keywords internal
#'
split_peak_names <- function(inp) {
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp),
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}



#' Make an ATAC cell dataset
#'
#' @keywords internal
#'
make_atac_cds <- function (input, binarize = FALSE) {

  if (is(input, "character")) {
    assertthat::is.readable(input)
    intersect_lean <- as.data.frame(data.table::fread(input,
                                                      header = FALSE))
  }
  else if (class(input) %in% c("matrix", "data.frame")) {
    intersect_lean <- input
  }
  else {
    stop("Input must be file path, matrix, or data.frame")
  }
  assertthat::assert_that(assertthat::are_equal(ncol(intersect_lean),
                                                3))
  assertthat::assert_that(is.logical(binarize))
  names(intersect_lean) <- c("site_name", "cell_name", "read_count")
  assertthat::assert_that(is.numeric(intersect_lean$read_count))
  intersect_lean$site_name <- as.factor(intersect_lean$site_name)
  intersect_lean$cell_name <- as.factor(intersect_lean$cell_name)
  cellinfo <- data.frame(cells = levels(intersect_lean$cell_name))
  row.names(cellinfo) <- cellinfo$cells
  cellinfo$temp <- seq_len(nrow(cellinfo))
  dhsinfo <- data.frame(site_name = levels(intersect_lean$site_name))
  dhsinfo <- cbind(dhsinfo, split_peak_names(dhsinfo$site_name))
  row.names(dhsinfo) <- dhsinfo$site_name
  names(dhsinfo) <- c("site_name", "chr", "bp1", "bp2")
  dhsinfo$chr <- gsub("chr", "", dhsinfo$chr)
  dhsinfo <- dhsinfo[order(as.character(dhsinfo$chr),
                           as.numeric(as.character(dhsinfo$bp2))),
  ]
  intersect_lean_ord <- intersect_lean[order(intersect_lean$site_name,
                                             intersect_lean$cell_name), ]
  dhsinfo <- dhsinfo[order(dhsinfo$site_name), ]
  cellinfo <- cellinfo[order(cellinfo$cells), ]
  intersect_lean_ord$site_name <- factor(intersect_lean_ord$site_name)
  intersect_lean_ord$cell_name <- factor(intersect_lean_ord$cell_name)
  intersect_lean_ord$site_name_num <- as.numeric(intersect_lean_ord$site_name)
  intersect_lean_ord$cell_name_num <- as.numeric(intersect_lean_ord$cell_name)
  if (binarize)
    intersect_lean_ord$read_count <- as.numeric(intersect_lean_ord$read_count >
                                                  0)
  sparse_intersect <- Matrix::sparseMatrix(i = intersect_lean_ord$site_name_num,
                                           j = intersect_lean_ord$cell_name_num,
                                           x = intersect_lean_ord$read_count)
  atac_cds <- suppressWarnings(monocle3::new_cell_data_set(methods::as(sparse_intersect,
                                                             "sparseMatrix"),
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = dhsinfo))
  pData(atac_cds)$temp <- NULL
  fData(atac_cds)$chr <- as.character(monocle3::fData(atac_cds)$chr)
  fData(atac_cds)$bp1 <- as.numeric(as.character(monocle3::fData(atac_cds)$bp1))
  fData(atac_cds)$bp2 <- as.numeric(as.character(monocle3::fData(atac_cds)$bp2))
  atac_cds <- atac_cds[order(monocle3::fData(atac_cds)$chr, monocle3::fData(atac_cds)$bp1),
  ]
  atac_cds <- monocle3::detect_genes(atac_cds)


  atac_cds
}



#' Run \code{cicero} to link coaccessibility enhancers
#'
#' @keywords internal
#'
run_cicero <- function (cds, genomic_coords, window = 5e+05, silent = FALSE,
                        sample_num = 100) {

  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(is.logical(silent))
  assertthat::assert_that(assertthat::is.number(window))
  assertthat::assert_that(assertthat::is.count(sample_num))
  if (!is.data.frame(genomic_coords)) {
    assertthat::is.readable(genomic_coords)
  }
  if (!silent)
    print("Starting Cicero")
  if (!silent)
    print("Calculating distance_parameter value")
  distance_parameters <- estimate_distance_parameter(cds,
                                                     window = window, maxit = 100,
                                                     sample_num = sample_num,
                                                     distance_constraint = 250000,
                                                     distance_parameter_convergence = 1e-22,
                                                     genomic_coords = genomic_coords)
  mean_distance_parameter <- mean(unlist(distance_parameters))
  if (!silent)
    print("Running models")
  cicero_out <- generate_cicero_models(cds, distance_parameter = mean_distance_parameter,
                                       window = window, genomic_coords = genomic_coords)
  if (!silent)
    print("Assembling connections")
  all_cons <- assemble_connections(cicero_out, silent = silent)
  if (!silent)
    print("Done")


  all_cons
}



#' Make a \code{cicero} cell dataset
#' @keywords internal
#'
make_cicero_cds <- function (cds, reduced_coordinates, k = 50, summary_stats = NULL,
                             size_factor_normalize = TRUE, silent = FALSE,
                             return_agg_info = FALSE) {

  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(is.data.frame(reduced_coordinates) |
                            is.matrix(reduced_coordinates))
  assertthat::assert_that(assertthat::are_equal(nrow(reduced_coordinates),
                                                nrow(monocle3::pData(cds))))
  assertthat::assert_that(setequal(row.names(reduced_coordinates),
                                   colnames(cds)))
  assertthat::assert_that(assertthat::is.count(k) & k > 1)
  assertthat::assert_that(is.character(summary_stats) | is.null(summary_stats))
  if (!is.null(summary_stats)) {
    assertthat::assert_that(all(summary_stats %in% names(monocle3::pData(cds))),
                            msg = paste("One of your summary_stats is missing",
                                        "from your pData table. Either add a", "column with the name in",
                                        "summary_stats, or remove the name", "from the summary_stats parameter.",
                                        collapse = " "))
    assertthat::assert_that(sum(vapply(summary_stats, function(x) {
      !(is(pData(cds)[, x], "numeric") | is(pData(cds)[,
                                                       x], "integer"))
    }, 1)) == 0, msg = paste("All columns in summary_stats must be",
                             "of class numeric or integer.", collapse = " "))
  }
  assertthat::assert_that(is.logical(size_factor_normalize))
  assertthat::assert_that(is.logical(silent))
  assertthat::assert_that(is.logical(return_agg_info))
  reduced_coordinates <- as.data.frame(reduced_coordinates)
  reduced_coordinates <- reduced_coordinates[colnames(cds),
  ]
  nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1,
                   replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0
  k2 <- k * 2
  get_shared <- function(other, this_choice) {
    k2 - length(union(cell_sample[other, ], this_choice))
  }
  while (length(good_choices) > 0 & it < 5000) {
    it <- it + 1
    choice <- sample(seq_len(length(good_choices)), size = 1,
                     replace = FALSE)
    new_chosen <- c(chosen, good_choices[choice])
    good_choices <- good_choices[good_choices != good_choices[choice]]
    cell_sample <- nn_map[new_chosen, ]
    others <- seq_len(nrow(cell_sample) - 1)
    this_choice <- cell_sample[nrow(cell_sample), ]
    shared <- sapply(others, get_shared, this_choice = this_choice)
    if (max(shared) < 0.9 * k) {
      chosen <- new_chosen
    }
  }
  cell_sample <- nn_map[chosen, ]
  if (!silent) {
    combs <- combn(nrow(cell_sample), 2)
    shared <- apply(combs, 2, function(x) {
      k2 - length(unique(as.vector(cell_sample[x, ])))
    })
    message(paste0("Overlap QC metrics:\nCells per bin: ",
                   k, "\nMaximum shared cells bin-bin: ", max(shared),
                   "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
                   median(shared)))
    if (mean(shared)/k > 0.1)
      warning("On average, more than 10% of cells are shared between paired bins.")
  }
  exprs_old <- monocle3::exprs(cds)
  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
                   cell_sample[x, , drop = FALSE])
  if (return_agg_info) {
    row.names(mask) <- colnames(exprs_old)
    colnames(mask) <- paste0("agg_", chosen)
    agg_map <- reshape2::melt(mask)
    agg_map <- agg_map[agg_map$value, ]
    agg_map$value <- NULL
    names(agg_map) <- c("cell", "agg_cell")
  }
  mask <- Matrix::Matrix(mask)
  new_exprs <- exprs_old %*% mask
  pdata <- monocle3::pData(cds)
  new_pcols <- "agg_cell"
  if (!is.null(summary_stats)) {
    new_pcols <- c(new_pcols, paste0("mean_", summary_stats))
  }
  new_pdata <- plyr::adply(cell_sample, 1, function(x) {
    sub <- pdata[x, ]
    df_l <- list()
    df_l["temp"] <- 1
    for (att in summary_stats) {
      df_l[paste0("mean_", att)] <- mean(sub[, att])
    }
    data.frame(df_l)
  })
  new_pdata$agg_cell <- paste("agg_", chosen, sep = "")
  new_pdata <- new_pdata[, new_pcols, drop = FALSE]
  row.names(new_pdata) <- new_pdata$agg_cell
  colnames(new_exprs) <- new_pdata$agg_cell
  fdf <- monocle3::fData(cds)
  new_pdata$temp <- NULL
  cicero_cds <- suppressWarnings(monocle3::new_cell_data_set(new_exprs,
                                                   cell_metadata = new_pdata, gene_metadata = fdf))
  cicero_cds <- monocle3::detect_genes(cicero_cds, min_expr = 0.1)
  cicero_cds <- estimate_size_factors(cicero_cds)
  if (any(!c("chr", "bp1", "bp2") %in% names(fData(cicero_cds)))) {
    fData(cicero_cds)$chr <- NULL
    fData(cicero_cds)$bp1 <- NULL
    fData(cicero_cds)$bp2 <- NULL
    fData(cicero_cds) <- cbind(monocle3::fData(cicero_cds),
                               df_for_coords(row.names(monocle3::fData(cicero_cds))))
  }
  if (size_factor_normalize) {
    cicero_cds <- suppressWarnings(monocle3::new_cell_data_set(Matrix::t(Matrix::t(monocle3::exprs(cicero_cds)) /
                                                                 monocle3::pData(cicero_cds)$Size_Factor),
                                                     cell_metadata = monocle3::pData(cicero_cds),
                                                     gene_metadata = monocle3::fData(cicero_cds)))
  }
  if (return_agg_info) {
    return(list(cicero_cds, agg_map))
  }


  cicero_cds
}



reconcile <- function(values) {
  if (length(values) == 1) return(values)
  if (sum(values >= 0) == length(values)) return(mean(values))
  if (sum(values <= 0) == length(values)) return(mean(values))
  if (sum(values == 0) == length(values)) return(0)
  return(NA_real_)
}

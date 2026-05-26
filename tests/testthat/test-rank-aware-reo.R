library(testthat)
library(LogicComm)

test_that("rank-aware evidence separates strength when binary LCS is saturated", {
  lr_db <- data.frame(
    lr_pair = "L_R", ligand = "L", receptor = "R",
    pathway = "rank_test", annotation = "synthetic",
    stringsAsFactors = FALSE
  )
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")

  make_rank <- function(value) {
    matrix(value, nrow = 2, ncol = 6,
           dimnames = list(c("L", "R"), paste0("C", seq_len(6))))
  }

  results <- list(
    Case1 = IdentifyRankLogicConsensus(rank_mat = make_rank(0.95), lr_db = lr_db,
                                       rank_threshold = 0.5,
                                       threshold_grid = c(0.5, 0.7, 0.9),
                                       verbose = FALSE),
    Case2 = IdentifyRankLogicConsensus(rank_mat = make_rank(0.90), lr_db = lr_db,
                                       rank_threshold = 0.5,
                                       threshold_grid = c(0.5, 0.7, 0.9),
                                       verbose = FALSE),
    Ctrl1 = IdentifyRankLogicConsensus(rank_mat = make_rank(0.60), lr_db = lr_db,
                                       rank_threshold = 0.5,
                                       threshold_grid = c(0.5, 0.7, 0.9),
                                       verbose = FALSE),
    Ctrl2 = IdentifyRankLogicConsensus(rank_mat = make_rank(0.58), lr_db = lr_db,
                                       rank_threshold = 0.5,
                                       threshold_grid = c(0.5, 0.7, 0.9),
                                       verbose = FALSE)
  )
  groups <- c(Case1 = "Case", Case2 = "Case", Ctrl1 = "Ctrl", Ctrl2 = "Ctrl")

  cmp <- CompareRankLogicGroups(results, group_info = groups, lr_db = lr_db,
                                verbose = FALSE)
  row <- cmp[cmp$lr_pair == "L_R", , drop = FALSE]

  expect_equal(results$Case1$binary_lcs, 1)
  expect_equal(results$Ctrl1$binary_lcs, 1)
  expect_gt(results$Case1$rank_dominance_lcs, results$Ctrl1$rank_dominance_lcs)
  expect_gt(row$delta_rank_score, 0.25)
  expect_equal(row$case_mean_binary_lcs, row$ctrl_mean_binary_lcs)
})

test_that("rank-aware specificity detects broad versus cell-type-pair concentrated axes", {
  lr_db <- data.frame(
    lr_pair = "L_R", ligand = "L", receptor = "R",
    pathway = "rank_test", annotation = "synthetic",
    stringsAsFactors = FALSE
  )
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")

  cells <- c("A1", "A2", "B1", "B2")
  labels <- c(A1 = "A", A2 = "A", B1 = "B", B2 = "B")
  knn <- Matrix::Matrix(0, nrow = 4, ncol = 4, sparse = TRUE)
  rownames(knn) <- colnames(knn) <- cells
  knn["A1", "B1"] <- 1
  knn["A2", "B2"] <- 1
  knn["B1", "A1"] <- 1
  knn["B2", "A2"] <- 1

  broad_rank <- matrix(0.9, nrow = 2, ncol = 4,
                       dimnames = list(c("L", "R"), cells))
  specific_rank <- matrix(0.1, nrow = 2, ncol = 4,
                          dimnames = list(c("L", "R"), cells))
  specific_rank["L", c("A1", "A2")] <- 0.9
  specific_rank["R", c("B1", "B2")] <- 0.9

  broad <- IdentifyRankLogicConsensus(rank_mat = broad_rank, knn_mat = knn,
                                      lr_db = lr_db, cell_labels = labels,
                                      rank_threshold = 0.5, verbose = FALSE)
  specific <- IdentifyRankLogicConsensus(rank_mat = specific_rank, knn_mat = knn,
                                         lr_db = lr_db, cell_labels = labels,
                                         rank_threshold = 0.5, verbose = FALSE)

  expect_lt(broad$specificity_score, specific$specificity_score)
  expect_equal(specific$top_specificity_group, "A|B")
  expect_equal(specific$n_specificity_groups, 1)
})

test_that("calc_REO_rank_score_matrix returns bounded aligned rank percentiles", {
  expr <- matrix(
    c(10, 1, 5,
      8, 2, 4,
      1, 9, 3),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("L", "R", "BG"), paste0("C", 1:3))
  )
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")

  rank_mat <- calc_REO_rank_score_matrix(expr, genes = c("L", "R"),
                                         verbose = FALSE)
  res <- IdentifyRankLogicConsensus(rank_mat = rank_mat, lr_db = lr_db,
                                    rank_threshold = 0.5, verbose = FALSE)

  expect_equal(rownames(rank_mat), c("L", "R"))
  expect_true(all(is.finite(rank_mat)))
  expect_true(all(rank_mat >= 0 & rank_mat <= 1))
  expect_true("rank_dominance_lcs" %in% names(res))
})


library(testthat)
library(LogicComm)

test_that("specificity and role interpretation helpers add publication fields", {
  reo <- Matrix::Matrix(
    matrix(c(
      1,1,1,0,0,0,
      1,1,0,1,1,0,
      0,0,1,1,1,1,
      1,0,1,0,1,0
    ), nrow = 4, byrow = TRUE,
    dimnames = list(c("L1", "R1", "L2", "R2"), paste0("C", 1:6))),
    sparse = TRUE
  )
  labels <- setNames(c("A", "A", "B", "B", "C", "C"), colnames(reo))
  knn <- Matrix::Matrix(1, nrow = 6, ncol = 6, sparse = TRUE)
  diag(knn) <- 0
  rownames(knn) <- colnames(knn) <- colnames(reo)
  lr_db <- data.frame(lr_pair = c("L1_R1", "L2_R2"), ligand = c("L1", "L2"),
                      receptor = c("R1", "R2"), pathway = c("TEST", "MHC-I"),
                      annotation = "unit", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L1", "L2")
  lr_db$receptor_genes <- list("R1", "R2")
  ct <- summarize_celltype_communication(reo, cell_labels = labels, knn_mat = knn,
                                         lr_db = lr_db, lcs_threshold = 0.01,
                                         min_edges = 1, verbose = FALSE)
  expect_true(all(c("role_separation_label", "communication_evidence_label",
                    "role_reliability_label", "dominant_role_strict",
                    "role_biological_interpretation") %in% names(ct$role_summary)))
  ct <- score_communication_specificity(ct, verbose = FALSE)
  expect_true("specificity_summary" %in% names(ct))
  expect_true(all(c("pair_specificity", "specificity_class", "identity_associated_flag") %in% names(ct$specificity_summary)))
  interp <- interpret_celltype_roles(ct)
  expect_true(is.data.frame(interp))
  expect_true("role_interpretation_caution" %in% names(interp))
})

test_that("permutation diagnostics handle degenerate positive null", {
  null_pair <- data.frame(sender_type = c("A", "B"), receiver_type = c("B", "A"),
                          observed = c(0.2, 0), null_mean = c(0, 0),
                          null_sd = c(0, 0.01), z_score = c(Inf, 0),
                          empirical_p = c(1/101, 1), fdr = c(0.02, 1),
                          n_perm = 100, n_null_nonmissing = 100,
                          min_possible_p = 1/101, min_possible_bh_fdr = 2/101,
                          degenerate_null = c(TRUE, FALSE),
                          degenerate_positive_null = c(TRUE, FALSE),
                          null_interpretation_flag = c("degenerate_positive_null", "low_resolution_null"),
                          metric = "sum_lcs", stringsAsFactors = FALSE)
  diag <- diagnose_permutation_resolution(null_pair)
  expect_true(inherits(diag, "LogicCommPermutationDiagnostic"))
  expect_equal(diag$n_degenerate_positive_null, 1)
})

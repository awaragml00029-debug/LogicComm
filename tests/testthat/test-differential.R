library(testthat)
library(LogicComm)

# Tests for subgroup-resolved multi-sample differential communication (v0.8.0).

make_diff_ct <- function(case) {
  reo <- Matrix::Matrix(
    matrix(c(if (case) 1 else 0, if (case) 1 else 0, 0, 0,  # L1 active in A only when case
             0, 0, 1, 1,                                    # R1 active in B
             1, 1, 0, 0,                                    # L2 active in A
             0, 0, 1, 1),                                   # R2 active in B
           nrow = 4, byrow = TRUE,
           dimnames = list(c("L1", "R1", "L2", "R2"), paste0("C", 1:4))),
    sparse = TRUE)
  labels <- stats::setNames(c("A", "A", "B", "B"), paste0("C", 1:4))
  lr_db <- data.frame(lr_pair = c("L1_R1", "L2_R2"), ligand = c("L1", "L2"),
                      receptor = c("R1", "R2"), pathway = c("PWa", "PWb"),
                      stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L1", "L2")
  lr_db$receptor_genes <- list("R1", "R2")
  summarize_celltype_communication(reo, cell_labels = labels, lr_db = lr_db, mode = "global",
                                   lcs_threshold = 0.1, min_edges = 1, verbose = FALSE)
}

make_diff <- function() {
  ct_list <- list(Case1 = make_diff_ct(TRUE), Case2 = make_diff_ct(TRUE),
                  Ctrl1 = make_diff_ct(FALSE), Ctrl2 = make_diff_ct(FALSE))
  gi <- c(Case1 = "Case", Case2 = "Case", Ctrl1 = "Ctrl", Ctrl2 = "Ctrl")
  differential_celltype_communication(ct_list, gi, lcs_threshold = 0.1,
                                      min_samples_per_group = 2, verbose = FALSE)
}

test_that("differential communication resolves the changed L-R pair to its sender->receiver subgroup", {
  diff <- make_diff()
  expect_s3_class(diff, "LogicCommDifferential")
  expect_true(all(c("sender_type", "receiver_type", "lr_pair", "asymmetry",
                    "direction", "pathway") %in% names(diff$lr)))

  ab <- diff$lr[diff$lr$sender_type == "A" & diff$lr$receiver_type == "B" &
                diff$lr$lr_pair == "L1_R1", ]
  expect_equal(nrow(ab), 1)
  expect_equal(ab$case_freq, 1)
  expect_equal(ab$ctrl_freq, 0)
  expect_equal(ab$asymmetry, 1)
  expect_equal(ab$direction, "Case-up")
  expect_equal(ab$pathway, "PWa")
})

test_that("subgroup summary points to the cell-type pair carrying the difference", {
  diff <- make_diff()
  sg <- diff$subgroup
  ab <- sg[sg$sender_type == "A" & sg$receiver_type == "B", ]
  expect_equal(ab$n_case_up, 1)
  expect_equal(ab$n_ctrl_up, 0)
  expect_equal(ab$dominant_direction, "Case-up")
  # No other subgroup should carry a Case/Ctrl difference in this design.
  expect_equal(sum(sg$n_changed), 1)
})

test_that("differential plots return ggplot objects and resolve subgroups", {
  diff <- make_diff()
  expect_s3_class(plot_differential_communication_summary(diff), "ggplot")
  expect_s3_class(plot_differential_celltype_heatmap(diff, metric = "asymmetry"), "ggplot")
  # The heatmap also accepts a raw CompareLogicGroups result with composite keys.
  expect_s3_class(plot_differential_celltype_heatmap(diff$comparison, metric = "asymmetry"), "ggplot")
})

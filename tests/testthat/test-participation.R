library(testthat)
library(LogicComm)

# Tests for the cell-type communication participation profile (v0.8.0).

make_participation_fixture <- function() {
  reo <- Matrix::Matrix(
    matrix(c(1, 1, 0, 0,   # L1 active in A cells
             0, 0, 1, 1,   # R1 active in B cells
             1, 1, 0, 0,   # L2 active in A cells
             0, 0, 1, 1),  # R2 active in B cells
           nrow = 4, byrow = TRUE,
           dimnames = list(c("L1", "R1", "L2", "R2"), paste0("C", 1:4))),
    sparse = TRUE
  )
  labels <- stats::setNames(c("A", "A", "B", "B"), paste0("C", 1:4))
  lr_db <- data.frame(
    lr_pair = c("L1_R1", "L2_R2"), ligand = c("L1", "L2"), receptor = c("R1", "R2"),
    pathway = c("PW_alpha", "PW_beta"), annotation = "t", stringsAsFactors = FALSE
  )
  lr_db$ligand_genes <- list("L1", "L2")
  lr_db$receptor_genes <- list("R1", "R2")
  ct <- summarize_celltype_communication(
    reo, cell_labels = labels, lr_db = lr_db, mode = "global",
    lcs_threshold = 0.1, min_edges = 1, verbose = FALSE
  )
  list(reo = reo, ct = ct)
}

test_that("summarize_celltype_participation reports sender/receiver fractions per cell type", {
  fx <- make_participation_fixture()
  part <- summarize_celltype_participation(fx$ct, fx$reo, verbose = FALSE)
  expect_s3_class(part, "LogicCommParticipation")

  gp <- part$group_participation
  a <- gp[gp$cell_type == "A", ]
  b <- gp[gp$cell_type == "B", ]
  # A cells express ligands only; B cells express receptors only.
  expect_equal(a$frac_sender_active, 1)
  expect_equal(a$frac_receiver_active, 0)
  expect_equal(a$frac_communicating, 1)
  expect_equal(b$frac_receiver_active, 1)
  expect_equal(b$frac_sender_active, 0)
  expect_true(all(gp$frac_communicating >= 0 & gp$frac_communicating <= 1))
  expect_true(all(c("importance_label", "importance_score", "total_active_lr_events") %in% names(gp)))
})

test_that("pathway and L-R composition shares sum to 1 within cell type and direction", {
  fx <- make_participation_fixture()
  part <- summarize_celltype_participation(fx$ct, fx$reo, verbose = FALSE)
  pc <- part$pathway_composition
  expect_true(all(c("cell_type", "direction", "pathway", "share") %in% names(pc)))
  agg <- tapply(pc$share, paste(pc$cell_type, pc$direction), sum)
  expect_true(all(abs(agg - 1) < 1e-8))

  lc <- part$lr_composition
  expect_true(all(c("lr_pair", "pathway", "share") %in% names(lc)))
})

test_that("participation plots return ggplot objects", {
  fx <- make_participation_fixture()
  part <- summarize_celltype_participation(fx$ct, fx$reo, verbose = FALSE)
  expect_s3_class(plot_celltype_participation(part), "ggplot")
  expect_s3_class(plot_celltype_participation(part, mode = "communicating"), "ggplot")
  expect_s3_class(plot_celltype_pathway_composition(part, direction = "outgoing"), "ggplot")
  expect_s3_class(plot_celltype_pathway_composition(part, direction = "incoming", level = "lr_pair"), "ggplot")
  expect_s3_class(plot_celltype_communication_profile(part, cell_type = "A"), "ggplot")
})

test_that("summarize_celltype_participation validates cell coverage", {
  fx <- make_participation_fixture()
  reo_missing <- fx$reo[, c("C1", "C2"), drop = FALSE]
  expect_error(
    summarize_celltype_participation(fx$ct, reo_missing, verbose = FALSE),
    "must contain all cells"
  )
})

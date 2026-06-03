library(testthat)
library(LogicComm)

# Tests for the cell-type-free gene-background ambient guard (v0.9.0).
# Reproduces the "ambient CD8A makes non-CD8 cells look CD8A+" artifact and
# checks that gene_background = "quantile" removes it without any cell labels.

make_ambient_fixture <- function(ambient = 3) {
  set.seed(11)
  nA <- 80; nB <- 80
  genes <- c("CD8A", "LIG", paste0("BG", 1:200))
  cells <- c(paste0("A", 1:nA), paste0("B", 1:nB))
  m <- matrix(rpois(length(genes) * length(cells), 0.4), length(genes), length(cells),
              dimnames = list(genes, cells))
  for (hk in paste0("BG", 1:6)) m[hk, ] <- rpois(ncol(m), 25)  # housekeeping (high)
  A <- 1:nA; B <- (nA + 1):(nA + nB)
  m["CD8A", A] <- rpois(nA, 10)       # type A genuinely expresses CD8A
  m["CD8A", B] <- rpois(nB, ambient)  # type B: ambient/contamination only
  m["LIG", ] <- rpois(ncol(m), 5)     # broadly expressed ligand
  list(m = m, labs = stats::setNames(c(rep("A", nA), rep("B", nB)), cells), A = A, B = B)
}

test_that("gene_background gate removes ambient receptor activity without cell-type labels", {
  d <- make_ambient_fixture(3)
  reo_none <- calc_REO_matrix(d$m, lr_genes = c("CD8A", "LIG"), gene_background = "none", verbose = FALSE)
  reo_gate <- calc_REO_matrix(d$m, lr_genes = c("CD8A", "LIG"), gene_background = "quantile", verbose = FALSE)

  b_active <- function(reo) mean(as.numeric(reo["CD8A", ])[d$B])
  a_active <- function(reo) mean(as.numeric(reo["CD8A", ])[d$A])

  # Ambient CD8A in type B is called active in many cells without the gate ...
  expect_gt(b_active(reo_none), 0.4)
  # ... and is largely removed with the gate.
  expect_lt(b_active(reo_gate), 0.15)
  # Genuine CD8A in type A is preserved.
  expect_gt(a_active(reo_gate), 0.6)

  expect_equal(attr(reo_gate, "gene_background"), "quantile")
  expect_equal(attr(reo_none, "gene_background"), "none")
})

test_that("gene_background gate flows into cell-type communication and removes the spurious axis", {
  d <- make_ambient_fixture(3)
  lr_db <- data.frame(lr_pair = "LIG_CD8A", ligand = "LIG", receptor = "CD8A",
                      pathway = "p", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("LIG")
  lr_db$receptor_genes <- list("CD8A")

  reo_none <- calc_REO_matrix(d$m, lr_genes = all_lr_genes(lr_db), gene_background = "none", verbose = FALSE)
  reo_gate <- calc_REO_matrix(d$m, lr_genes = all_lr_genes(lr_db), gene_background = "quantile", verbose = FALSE)
  ct_none <- summarize_celltype_communication(reo_none, cell_labels = d$labs, lr_db = lr_db,
                                              mode = "global", lcs_threshold = 0.01, min_edges = 1, verbose = FALSE)
  ct_gate <- summarize_celltype_communication(reo_gate, cell_labels = d$labs, lr_db = lr_db,
                                              mode = "global", lcs_threshold = 0.01, min_edges = 1, verbose = FALSE)
  ab_none <- ct_none$lr_table$lcs[ct_none$lr_table$sender_type == "A" & ct_none$lr_table$receiver_type == "B"]
  ab_gate <- ct_gate$lr_table$lcs[ct_gate$lr_table$sender_type == "A" & ct_gate$lr_table$receiver_type == "B"]

  expect_gt(ab_none, 0.4)   # spurious A -> B via ambient CD8A is high by default
  expect_lt(ab_gate, 0.1)   # and is gated out, with no cell-type labels involved in the gate
})

test_that("gene_background validates its quantile argument", {
  d <- make_ambient_fixture(1)
  expect_error(calc_REO_matrix(d$m, lr_genes = "CD8A", gene_background = "quantile",
                               gene_background_quantile = 1.5, verbose = FALSE))
  expect_error(calc_REO_matrix(d$m, lr_genes = "CD8A", gene_background = "bogus", verbose = FALSE))
})

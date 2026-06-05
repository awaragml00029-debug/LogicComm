library(testthat)
library(LogicComm)

# Regression tests for diagnose_proliferation_confound() (v0.9.3). A cycling cell
# type has both a high cell-cycle signature and a broad transcriptome (many
# detected genes); the diagnostic should flag it as a proliferation/breadth hub
# and annotate axes touching it, without altering any scoring columns.

make_prolif_fixture <- function(seed = 1) {
  set.seed(seed)
  types <- c("Cycling", paste0("Rest", 1:5))
  per <- 40L
  cells <- paste0("c", seq_len(per * length(types)))
  labels <- stats::setNames(rep(types, each = per), cells)
  cyc_genes <- paste0("CC", 1:12)
  other <- paste0("G", 1:300)
  genes <- c(cyc_genes, other)
  m <- matrix(0L, length(genes), length(cells), dimnames = list(genes, cells))
  for (j in seq_along(labels)) {
    if (labels[j] == "Cycling") {
      on <- sample(other, 180)
      m[on, j] <- rpois(length(on), 3)
      m[cyc_genes, j] <- rpois(length(cyc_genes), 15)    # strong cycle signal + broad
    } else {
      on <- sample(other, 45)
      m[on, j] <- rpois(length(on), 2)
      m[cyc_genes, j] <- rpois(length(cyc_genes), 0.05)  # ambient cycle genes only
    }
  }
  lr <- data.frame(
    sender_type   = c("Cycling", "Rest1"),
    receiver_type = c("Rest1",   "Rest2"),
    lr_pair = c("L1_R1", "L2_R2"),
    active = c(TRUE, TRUE), lcs = c(0.5, 0.5),
    stringsAsFactors = FALSE
  )
  ct <- list(cell_labels = labels, lr_table = lr,
             role_summary = data.frame(cell_type = types, stringsAsFactors = FALSE))
  class(ct) <- "LogicCommCellTypeComm"
  list(ct = ct, expr = m, cyc_genes = cyc_genes)
}

test_that("the cycling cell type is flagged as a proliferation/breadth hub", {
  f <- make_prolif_fixture()
  out <- diagnose_proliferation_confound(
    f$ct, f$expr, s_genes = character(0), g2m_genes = f$cyc_genes, verbose = FALSE)

  ps <- out$proliferation_summary
  hub <- ps$cell_type[ps$proliferation_hub_flag %in% TRUE]
  expect_identical(hub, "Cycling")

  cy <- ps[ps$cell_type == "Cycling", ]
  re <- ps[ps$cell_type == "Rest1", ]
  expect_gt(cy$mean_cycling_score, re$mean_cycling_score)
  expect_gt(cy$mean_breadth, re$mean_breadth)
})

test_that("axes touching the hub are annotated; scoring is untouched", {
  f <- make_prolif_fixture()
  out <- diagnose_proliferation_confound(
    f$ct, f$expr, s_genes = character(0), g2m_genes = f$cyc_genes, verbose = FALSE)
  lr <- out$lr_table

  sender_hub <- lr[lr$sender_type == "Cycling", ]
  clean <- lr[lr$sender_type == "Rest1" & lr$receiver_type == "Rest2", ]
  expect_true(sender_hub$proliferation_confound_flag)
  expect_identical(sender_hub$proliferation_hub_role, "sender")
  expect_false(clean$proliferation_confound_flag)
  expect_identical(clean$proliferation_hub_role, "none")

  # Diagnostic only: it must not introduce scoring columns.
  expect_false("discovery_score" %in% names(lr))
  expect_false("evidence_tier" %in% names(lr))
})

test_that("expr / cell-label mismatch is an error", {
  f <- make_prolif_fixture()
  bad <- f$expr[, 1:10]
  colnames(bad) <- NULL
  expect_error(
    diagnose_proliferation_confound(f$ct, bad, g2m_genes = f$cyc_genes, verbose = FALSE)
  )
})

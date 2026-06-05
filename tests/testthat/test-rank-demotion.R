library(testthat)
library(LogicComm)

# Regression tests for broad/ubiquitous axis demotion in rank_communication_axes()
# (v0.9.3). score_communication_specificity() flags broad axes, but until 0.9.3
# the flag was annotation only and a strong, null-supported broad axis still
# reached Tier 1/2. demote_broad caps such axes at Tier 3 and penalizes their
# discovery_score, without touching the active flag or LCS.

make_rank_fixture <- function() {
  lr <- data.frame(
    sender_type = c("A", "B", "C"),
    receiver_type = c("X", "Y", "Z"),
    lr_pair = c("HLA-A_CD8A", "FOO_BAR", "BAZ_QUX"),
    ligand = c("HLA-A", "FOO", "BAZ"),
    receptor = c("CD8A", "BAR", "QUX"),
    pathway = c("MHC-I", "P2", "P3"),
    lcs = c(0.9, 0.8, 0.2),
    active = c(TRUE, TRUE, TRUE),
    pair_specificity = c(0.10, 0.90, 0.50),
    ubiquitous_interaction_flag = c(TRUE, FALSE, FALSE),
    identity_associated_flag = c(TRUE, FALSE, FALSE),
    specificity_class = c("Identity-associated broad axis",
                          "Pair-specific candidate", "Intermediate specificity"),
    stringsAsFactors = FALSE
  )
  ct <- list(lr_table = lr,
             specificity_summary = data.frame(lr_pair = lr$lr_pair,
                                              stringsAsFactors = FALSE))
  class(ct) <- "LogicCommCellTypeComm"
  ct
}

get_row <- function(d, s) d[d$sender_type == s, , drop = FALSE]

test_that("demote_broad caps broad axes at Tier 3 and penalizes discovery_score", {
  ct <- make_rank_fixture()
  np <- data.frame(sender_type = c("A", "B", "C"),
                   receiver_type = c("X", "Y", "Z"),
                   empirical_p = c(0, 0, 0), stringsAsFactors = FALSE)

  on  <- rank_communication_axes(ct, null_pair = np)                  # default TRUE
  off <- rank_communication_axes(ct, null_pair = np, demote_broad = FALSE)

  # The broad axis is detected and flagged.
  expect_true(get_row(on, "A")$broad_axis_flag)
  expect_false(get_row(on, "B")$broad_axis_flag)

  # Demotion ON: broad axis capped at Tier 3 with a dedicated label; the genuine
  # specific axis keeps its high tier.
  expect_identical(get_row(on, "A")$evidence_tier, "Tier 3: broad / non-specific")
  expect_identical(get_row(on, "B")$evidence_tier, "Tier 1: strong, specific, supported")

  # discovery_score of the broad axis is exactly broad_penalty x the un-demoted score.
  expect_equal(get_row(on, "A")$discovery_score,
               0.5 * get_row(off, "A")$discovery_score)

  # Demotion OFF: the broad axis would otherwise be a Tier 2 "strong candidate".
  expect_identical(get_row(off, "A")$evidence_tier, "Tier 2: strong candidate")

  # The specific axis outranks the demoted broad axis.
  expect_identical(on$sender_type[1], "B")
  expect_gt(which(on$sender_type == "A"), which(on$sender_type == "B"))
})

test_that("demote_broad does not change the active flag or lcs", {
  ct <- make_rank_fixture()
  on <- rank_communication_axes(ct)
  expect_true(all(on$active))
  expect_setequal(on$lcs, c(0.9, 0.8, 0.2))
})

test_that("broad_penalty is validated", {
  ct <- make_rank_fixture()
  expect_error(rank_communication_axes(ct, broad_penalty = 1.5))
  expect_error(rank_communication_axes(ct, broad_penalty = -0.1))
})

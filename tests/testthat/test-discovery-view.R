library(testthat)
library(LogicComm)

# Tests for communication_discovery_view() (v0.9.4). It filters a ranked table to
# the confound-free candidates and, by default, re-scales strength/specificity
# and re-tiers WITHIN the kept set so genuine but lower-LCS axes are tiered on
# their own merits instead of against the broad/cycling noise that was removed.

make_view_fixture <- function() {
  data.frame(
    sender_type   = c("A", "Cyc", "NK", "CD8", "Treg"),
    receiver_type = c("B", "X", "NK", "Treg", "NK"),
    lr_pair       = c("HLA-B_CD8A", "FOO_BAR", "HLA-E_KLRC1", "TNFSF9_TNFRSF9", "CCL5_CCR4"),
    lcs                         = c(0.90, 0.80, 0.30, 0.40, 0.05),
    pair_specificity            = c(0.10, 0.20, 0.50, 0.95, 0.20),
    threshold_stability         = c(1.0, 1.0, 0.8, 1.0, 0.3),
    null_support                = c(0.9, 0.9, 0.7, 0.95, 0.2),
    broad_axis_flag             = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    proliferation_confound_flag = c(FALSE, TRUE, FALSE, FALSE, FALSE),
    identity_associated_flag    = c(TRUE, FALSE, TRUE, FALSE, FALSE),
    discovery_score             = c(0.70, 0.65, 0.55, 0.50, 0.30),
    evidence_tier = c("Tier 3: broad / non-specific", "Tier 2: strong candidate",
                      "Tier 3: emerging candidate", "Tier 3: emerging candidate",
                      "Tier 4: weak / context-dependent"),
    stringsAsFactors = FALSE
  )
}

test_that("default view drops broad + proliferation axes and re-tiers the survivors", {
  v <- communication_discovery_view(make_view_fixture(), verbose = FALSE)

  # Broad (HLA-B_CD8A) and proliferation (FOO_BAR) axes are removed.
  expect_equal(nrow(v), 3)
  expect_false(any(v$lr_pair %in% c("HLA-B_CD8A", "FOO_BAR")))

  # TNFSF9_TNFRSF9 is the strongest + most specific + supported survivor, so after
  # re-scaling within the kept set it is promoted to Tier 1 and ranks first --
  # even though it was only Tier 3 against the full noisy set.
  expect_equal(v$lr_pair[1], "TNFSF9_TNFRSF9")
  expect_equal(v$evidence_tier[v$lr_pair == "TNFSF9_TNFRSF9"],
               "Tier 1: strong, specific, supported")
  expect_true(all(diff(v$discovery_score) <= 0))
})

test_that("drop_identity also removes identity-associated axes", {
  v <- communication_discovery_view(make_view_fixture(), drop_identity = TRUE,
                                    verbose = FALSE)
  expect_false("HLA-E_KLRC1" %in% v$lr_pair)
  expect_equal(nrow(v), 2)                       # TNFSF9 + CCL5 remain
})

test_that("rescale = FALSE keeps the original scores and only filters + orders", {
  fx <- make_view_fixture()
  v <- communication_discovery_view(fx, rescale = FALSE, verbose = FALSE)
  expect_equal(nrow(v), 3)
  # Original discovery_score preserved (not recomputed).
  expect_equal(v$discovery_score[v$lr_pair == "HLA-E_KLRC1"], 0.55)
  expect_equal(v$evidence_tier[v$lr_pair == "HLA-E_KLRC1"], "Tier 3: emerging candidate")
  expect_true(all(diff(v$discovery_score) <= 0))
})

test_that("missing proliferation column warns but still drops broad axes", {
  fx <- make_view_fixture()
  fx$proliferation_confound_flag <- NULL
  expect_warning(
    v <- communication_discovery_view(fx, verbose = FALSE),
    "proliferation_confound_flag"
  )
  expect_false("HLA-B_CD8A" %in% v$lr_pair)      # broad still removed
})

test_that("top_n limits the returned rows and a non-ranked input errors", {
  v <- communication_discovery_view(make_view_fixture(), top_n = 1, verbose = FALSE)
  expect_equal(nrow(v), 1)
  expect_equal(v$lr_pair, "TNFSF9_TNFRSF9")
  expect_error(communication_discovery_view(data.frame(x = 1)), "rank_communication_axes")
})

test_that("rank_communication_axes is unchanged by the shared-core refactor", {
  # A minimal hand-built ct_comm exercising the demotion path through the ranker.
  lr <- data.frame(
    sender_type = c("A", "B", "C"), receiver_type = c("X", "Y", "Z"),
    lr_pair = c("HLA-A_CD8A", "FOO_BAR", "BAZ_QUX"),
    lcs = c(0.9, 0.8, 0.2), active = c(TRUE, TRUE, TRUE),
    pair_specificity = c(0.1, 0.9, 0.5),
    ubiquitous_interaction_flag = c(TRUE, FALSE, FALSE),
    identity_associated_flag = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  ct <- list(lr_table = lr, specificity_summary = data.frame(x = 1))
  class(ct) <- "LogicCommCellTypeComm"
  null_pair <- data.frame(sender_type = c("A", "B", "C"),
                          receiver_type = c("X", "Y", "Z"),
                          empirical_p = c(0, 0, 0))
  r <- rank_communication_axes(ct, null_pair = null_pair)

  # Broad axis is demoted, never Tier 1/2; specific axis tops the ranking.
  expect_true("broad_axis_flag" %in% names(r))
  expect_equal(r$evidence_tier[r$lr_pair == "HLA-A_CD8A"], "Tier 3: broad / non-specific")
  expect_false(any(r$broad_axis_flag & grepl("Tier 1|Tier 2", r$evidence_tier)))
  expect_equal(r$lr_pair[1], "FOO_BAR")
})

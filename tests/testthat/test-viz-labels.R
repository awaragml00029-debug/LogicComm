library(testthat)
library(LogicComm)

# Tests for publication-figure label readability fixes (v0.8.0).

test_that("label helpers shorten long names", {
  expect_equal(.short_label("short", 30L), "short")
  long <- paste(rep("X", 50), collapse = "")
  out <- .short_label(long, 20L)
  expect_true(nchar(out) <= 20L)
  expect_match(out, "\u2026$")

  key <- "Activated_Cytotoxic_CD8_T_GZMK_NKG7|Tissue_NK_XCL1_XCL2_AREG|VEGFA_KDR"
  cl <- .compact_feature_label(key, 40L)
  expect_true(startsWith(cl, "VEGFA_KDR"))   # L-R pair leads the label
  expect_true(nchar(cl) <= 40L)
  expect_match(cl, "\u2192")                 # sender -> receiver arrow

  # A non-composite value is returned (truncated) unchanged in structure.
  expect_equal(.compact_feature_label("VEGFA_KDR", 40L), "VEGFA_KDR")
})

test_that("plot_specificity_stability honours label_field and shortens labels", {
  sens <- list(stability = data.frame(
    feature = c("SenderTypeWithAVeryLongName|ReceiverTypeAlsoVeryLong|VEGFA_KDR",
                "SenderTypeWithAVeryLongName|ReceiverTypeAlsoVeryLong|TGFB1_TGFBR1"),
    active_fraction = c(1, 0.5), stringsAsFactors = FALSE))
  specificity <- data.frame(
    lr_pair = c("VEGFA_KDR", "TGFB1_TGFBR1"),
    pair_specificity = c(0.9, 0.3),
    pathway = c("VEGF", "TGFb"),
    specificity_class = c("Pair-specific candidate", "Broad/ubiquitous axis"),
    total_lcs = c(2, 1), stringsAsFactors = FALSE)

  p_lr <- plot_specificity_stability(sens, specificity = specificity, label_field = "lr_pair")
  expect_s3_class(p_lr, "ggplot")
  expect_s3_class(ggplot2::ggplotGrob(p_lr), "gtable")   # renders without error
  lab_lr <- p_lr$layers[[2]]$data$label
  expect_true("VEGFA_KDR" %in% lab_lr)                    # short L-R name, not the full key

  p_feat <- plot_specificity_stability(sens, specificity = specificity, label_field = "feature")
  expect_s3_class(ggplot2::ggplotGrob(p_feat), "gtable")
  expect_true(all(nchar(p_feat$layers[[2]]$data$label) <= 40L))  # compact, not the raw key

  p_none <- plot_specificity_stability(sens, specificity = specificity, label_field = "none")
  expect_s3_class(p_none, "ggplot")
  expect_false(any(p_none$data$to_label))                # nothing labelled
})

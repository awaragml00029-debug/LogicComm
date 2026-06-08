test_that("calc_rank_shift uses the same analyzable samples for ranks and expression means", {
  make_sample <- function(g_values) {
    matrix(
      c(g_values,
        rep(5, length(g_values))),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("G", "BG"), paste0("C", seq_along(g_values)))
    )
  }

  sample_list <- list(
    Case1 = make_sample(c(10, 10, 10, 10)),
    Ctrl1 = make_sample(c(1, 1, 1, 1)),
    CtrlFiltered = make_sample(c(100, 0, 0, 0))
  )
  group_info <- c(Case1 = "Case", Ctrl1 = "Ctrl", CtrlFiltered = "Ctrl")
  lr_db <- data.frame(lr_pair = "G_R", ligand = "G", receptor = "R",
                      stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("G")
  lr_db$receptor_genes <- list("R")

  result <- calc_rank_shift(
    sample_list,
    group_info = group_info,
    genes = "G",
    lr_db = lr_db,
    min_detection_frac = 0.5,
    verbose = FALSE
  )

  row <- result[result$gene == "G", , drop = FALSE]
  expect_equal(row$mean_expr_case, 10)
  expect_equal(row$mean_expr_ctrl, 1)
  expect_equal(row$log2fc_expr, round(log2(10.1 / 1.1), 3))
})

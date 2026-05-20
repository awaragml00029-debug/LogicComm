test_that("logic_* core API prepares, scores, and summarizes communication", {
  counts <- Matrix::Matrix(
    matrix(c(10, 0, 10,
             0, 10, 10,
             1, 1, 1), nrow = 3, byrow = TRUE,
           dimnames = list(c("L", "R", "BG"), c("C1", "C2", "C3"))),
    sparse = TRUE
  )
  lr_db <- data.frame(
    lr_pair = "L_R",
    ligand = "L",
    receptor = "R",
    pathway = "test",
    annotation = "test",
    stringsAsFactors = FALSE
  )
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")

  check <- logic_check_lrdb(lr_db)
  expect_s3_class(check, "LogicCommLRDBCheck")
  expect_true(check$ok)
  expect_equal(logic_get_lr_genes(lr_db), c("L", "R"))

  reo <- logic_prepare(counts, lr_db = lr_db, rank_threshold = 0.5, verbose = FALSE)
  expect_true(inherits(reo, "sparseMatrix"))
  expect_equal(rownames(reo), c("L", "R"))

  reo_lcs <- Matrix::Matrix(
    matrix(c(1, 0, 1,
             0, 1, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), c("C1", "C2", "C3"))),
    sparse = TRUE
  )
  lcs <- logic_score_lr(reo_lcs, lr_db = lr_db, mode = "global", verbose = FALSE)
  expect_s3_class(lcs, "LCSVector")
  expect_equal(unname(lcs["L_R"]), 1 / 3)

  labels <- setNames(c("A", "B", "B"), colnames(reo_lcs))
  ct <- logic_summarize_celltypes(
    reo_lcs,
    cell_labels = labels,
    lr_db = lr_db,
    mode = "global",
    lcs_threshold = 0.01,
    min_edges = 1,
    verbose = FALSE
  )
  expect_s3_class(ct, "LogicCommCellTypeComm")
  expect_true(all(c("lr_table", "pair_summary", "role_summary") %in% names(ct)))
})

test_that("logic_compare_groups and logic_run expose clear group workflow", {
  lcs_list <- list(
    S1 = setNames(c(0.2, 0.0), c("L_R", "X_Y")),
    S2 = setNames(c(0.1, 0.0), c("L_R", "X_Y")),
    S3 = setNames(c(0.0, 0.3), c("L_R", "X_Y")),
    S4 = setNames(c(0.0, 0.2), c("L_R", "X_Y"))
  )
  groups <- c(S1 = "Case", S2 = "Case", S3 = "Ctrl", S4 = "Ctrl")

  cmp <- logic_compare_groups(lcs_list, group_info = groups, verbose = FALSE)
  expect_s3_class(cmp, "LogicCommResult")
  expect_true("asymmetry" %in% names(cmp))

  counts1 <- Matrix::Matrix(
    matrix(c(10, 0, 10,
             0, 10, 10,
             1, 1, 1), nrow = 3, byrow = TRUE,
           dimnames = list(c("L", "R", "BG"), c("C1", "C2", "C3"))),
    sparse = TRUE
  )
  counts2 <- counts1
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  samples <- list(S1 = counts1, S2 = counts2)
  groups2 <- c(S1 = "Case", S2 = "Ctrl")

  analysis <- logic_run(samples, group_info = groups2, lr_db = lr_db,
                        rank_threshold = 0.5, verbose = FALSE)
  expect_s3_class(analysis, "LogicCommAnalysis")
  expect_s3_class(analysis, "LogicCommMulti")
  expect_s3_class(logic_compare_groups(analysis), "LogicCommResult")
})

test_that("logic_score_lr validates requested neighborhood mode", {
  reo <- Matrix::Matrix(
    matrix(c(1, 0,
             0, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), c("C1", "C2"))),
    sparse = TRUE
  )
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")

  expect_error(
    logic_score_lr(reo, lr_db = lr_db, mode = "neighborhood", verbose = FALSE),
    "requires seurat_obj or knn_mat"
  )
})

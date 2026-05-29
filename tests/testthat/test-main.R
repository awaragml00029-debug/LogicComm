library(testthat)
library(LogicComm)

# ── Generate synthetic data ───────────────────────────────────────────────────
set.seed(42)
n_genes <- 50
n_cells <- 200

gene_names <- c(
  "VEGFA","KDR","FLT1","TGFB1","TGFBR1","TGFBR2",
  "CXCL12","CXCR4","TNF","TNFRSF1A","IL6","IL6R",
  "COL1A1","CD44","PDGFB","PDGFRB","FGF2","FGFR1",
  paste0("GENE", seq_len(n_genes - 18))
)
cell_names <- paste0("Cell", seq_len(n_cells))

# Simulate sparse count matrix
count_mat <- matrix(
  stats::rpois(n_genes * n_cells, lambda = 1.5),
  nrow = n_genes, ncol = n_cells,
  dimnames = list(gene_names, cell_names)
)
# Make some genes highly expressed in first 100 cells (simulated signal)
count_mat["VEGFA", 1:100]  <- count_mat["VEGFA", 1:100] + 5
count_mat["KDR",   1:100]  <- count_mat["KDR",   1:100] + 4


# ── Tests ─────────────────────────────────────────────────────────────────────

test_that("calc_REO_matrix works on numeric matrix", {
  
  lr_g <- all_lr_genes(lr_pairs_human)

  reo <- calc_REO_matrix(count_mat, lr_genes = lr_g, verbose = FALSE)
  expect_true(inherits(reo, "sparseMatrix"))
  expect_true(all(reo@x %in% c(0L, 1L)))
  expect_equal(ncol(reo), n_cells)
})

test_that("calc_REO_matrix handles chunk processing", {
  
  lr_g <- all_lr_genes(lr_pairs_human)

  reo1 <- calc_REO_matrix(count_mat, lr_genes = lr_g, chunk_size = 50,  verbose = FALSE)
  reo2 <- calc_REO_matrix(count_mat, lr_genes = lr_g, chunk_size = 500, verbose = FALSE)
  expect_equal(as.matrix(reo1), as.matrix(reo2))
})

test_that("IdentifyLogicConsensus global mode returns named numeric vector", {
  
  lr_g <- all_lr_genes(lr_pairs_human)
  reo  <- calc_REO_matrix(count_mat, lr_genes = lr_g, verbose = FALSE)

  lcs  <- IdentifyLogicConsensus(reo, verbose = FALSE)
  expect_true(is.numeric(lcs))
  expect_true(all(lcs[!is.na(lcs)] >= 0 & lcs[!is.na(lcs)] <= 1))
  expect_true("VEGFA_KDR" %in% names(lcs))
})

test_that("CompareLogicGroups returns LogicCommResult", {
  
  lr_g <- all_lr_genes(lr_pairs_human)

  # Simulate 6 samples: 3 Case + 3 Ctrl
  make_lcs <- function(seed, boost = 0) {
    set.seed(seed)
    mat_i <- count_mat
    mat_i["VEGFA", ] <- mat_i["VEGFA", ] + boost
    reo_i <- calc_REO_matrix(mat_i, lr_genes = lr_g, verbose = FALSE)
    IdentifyLogicConsensus(reo_i, verbose = FALSE)
  }

  lcs_list <- c(
    lapply(1:3, function(s) make_lcs(s, boost = 5)),  # Case: boosted VEGFA
    lapply(4:6, function(s) make_lcs(s, boost = 0))   # Ctrl: baseline
  )
  names(lcs_list) <- c(paste0("Case_", 1:3), paste0("Ctrl_", 1:3))
  groups <- c(rep("Case", 3), rep("Ctrl", 3))
  names(groups) <- names(lcs_list)

  result <- CompareLogicGroups(lcs_list, group_info = groups, verbose = FALSE)

  expect_true(inherits(result, "LogicCommResult"))
  expect_true(all(c("lr_pair","case_freq","ctrl_freq","asymmetry","fdr_fisher") %in% names(result)))
  expect_true(nrow(result) > 0)
  # VEGFA-KDR should be enriched in Case
  vegfa_row <- result[result$lr_pair == "VEGFA_KDR", ]
  expect_true(nrow(vegfa_row) == 1)
  expect_true(vegfa_row$case_freq >= vegfa_row$ctrl_freq)
})

test_that("run_multisample validates four named PBMC-style samples with a known perturbation", {
  genes <- c("MIF", "CD74", paste0("BG", seq_len(20)))
  make_sample <- function(active = FALSE) {
    mat <- matrix(1, nrow = length(genes), ncol = 12,
                  dimnames = list(genes, paste0("C", seq_len(12))))
    if (isTRUE(active)) {
      mat["MIF", ] <- 20
      mat["CD74", ] <- 20
    } else {
      mat["MIF", ] <- 0
      mat["CD74", ] <- 0
    }
    Matrix::Matrix(mat, sparse = TRUE)
  }

  lr_db <- data.frame(lr_pair = "MIF_CD74", ligand = "MIF", receptor = "CD74",
                      pathway = "MIF", annotation = "demo", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("MIF")
  lr_db$receptor_genes <- list("CD74")

  samples <- list(
    pbmc1 = make_sample(TRUE),
    pbmc2 = make_sample(TRUE),
    pbmc3 = make_sample(FALSE),
    pbmc4 = make_sample(FALSE)
  )
  groups <- c(pbmc1 = "Case", pbmc2 = "Case", pbmc3 = "Ctrl", pbmc4 = "Ctrl")

  result <- run_multisample(
    samples,
    group_info = groups,
    lr_db = lr_db,
    rank_threshold = 0.5,
    lcs_threshold = 0.5,
    min_samples_per_group = 2,
    verbose = FALSE
  )

  row <- result$comparison[result$comparison$lr_pair == "MIF_CD74", , drop = FALSE]
  expect_s3_class(result, "LogicCommMulti")
  expect_equal(names(result$lcs_list), names(samples))
  expect_equal(row$case_freq, 1)
  expect_equal(row$ctrl_freq, 0)
  expect_equal(row$asymmetry, 1)
  expect_gt(row$case_mean_lcs, row$ctrl_mean_lcs)
})

test_that("score_lr_activity returns valid scores", {
  
  lr_g <- all_lr_genes(lr_pairs_human)
  reo  <- calc_REO_matrix(count_mat, lr_genes = lr_g, verbose = FALSE)

  scores <- score_lr_activity(reo, verbose = FALSE)
  expect_true(is.list(scores))
  expect_true(all(scores$comm_score >= 0 & scores$comm_score <= 1))
  expect_equal(length(scores$comm_score), n_cells)
})

test_that("AutoLabelLogicClusters returns data frame", {
  
  lr_g    <- all_lr_genes(lr_pairs_human)
  reo     <- calc_REO_matrix(count_mat, lr_genes = lr_g, verbose = FALSE)
  clusters <- setNames(rep(c("TypeA","TypeB"), each = n_cells/2), cell_names)

  labels <- AutoLabelLogicClusters(
    reo, cell_labels = clusters,
    lr_pairs = c("VEGFA_KDR","TGFB1_TGFBR1"),
    verbose = FALSE)

  expect_true(is.data.frame(labels))
  expect_true("sender_cluster" %in% names(labels))
})

test_that("all_lr_genes returns character vector", {

  genes <- all_lr_genes(lr_pairs_human)
  expect_true(is.character(genes))
  expect_true(length(genes) > 10)
  expect_true("VEGFA" %in% genes)
})

test_that("CellChatDB conversion resolves simple, complex, and modulator interactions", {
  cellchat_db <- list(
    interaction = data.frame(
      interaction_name = c("TGFB1_TGFBR1_TGFBR2", "CXCL12_CXCR4"),
      ligand = c("TGFB1", "CXCL12"),
      receptor = c("TGFbR_complex", "CXCR4"),
      pathway_name = c("TGFb", "CXCL"),
      agonist = c("TGFBR3", ""),
      antagonist = c("", "ACKR3"),
      co_A_receptor = c("", ""),
      co_I_receptor = c("", ""),
      stringsAsFactors = FALSE
    ),
    complex = data.frame(
      subunit_1 = "TGFBR1",
      subunit_2 = "TGFBR2",
      stringsAsFactors = FALSE,
      row.names = "TGFbR_complex"
    )
  )

  converted <- as_logiccomm_lr_db_from_cellchat(cellchat_db)
  expect_true(all(c("lr_pair", "ligand_genes", "receptor_genes", "pathway") %in% names(converted)))
  expect_true(all(c("agonist_genes", "antagonist_genes", "co_A_receptor_genes", "co_I_receptor_genes") %in% names(converted)))
  expect_equal(converted$ligand_genes[[1]], "TGFB1")
  expect_equal(converted$receptor_genes[[1]], c("TGFBR1", "TGFBR2"))
  expect_equal(converted$agonist_genes[[1]], "TGFBR3")
  expect_equal(converted$antagonist_genes[[2]], "ACKR3")
  expect_equal(converted$annotation[[1]], "CellChatDB")
})

test_that("rank-aware logic gate blocks antagonist rank reversal", {
  expr <- matrix(
    c(10, 10, 10, 10,
      8,  8,  8,  8,
      1,  1,  1,  1,
      1,  9,  1,  9),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("L", "R", "LOW", "ANT"), paste0("C", 1:4))
  )
  reo_res <- calc_REO_matrix(expr, lr_genes = c("L", "R", "ANT"),
                             rank_threshold = 0.25, return_rank = TRUE,
                             verbose = FALSE)
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
  lr_db$antagonist_genes <- list("ANT")
  lr_db$agonist_genes <- list(character(0))
  lr_db$co_A_receptor_genes <- list(character(0))
  lr_db$co_I_receptor_genes <- list(character(0))

  gate <- IdentifyLogicGateConsensus(reo_res, lr_db = lr_db,
                                     negative_gate = "rank_block",
                                     antagonist_side = "receiver",
                                     verbose = FALSE)
  expect_equal(gate$base_lcs, 1)
  expect_equal(gate$regulated_lcs, 0.5)
  expect_equal(gate$block_rate, 0.5)

  summary <- summarize_lr_modulators(reo_res$logic, lr_db, rank_mat = reo_res$rank)
  expect_true(summary$has_antagonist)
  expect_equal(summary$antagonist_outranks_receptor_rate, 0.5)
})

test_that("KNN mode validates and aligns cell names", {
  reo <- Matrix::Matrix(
    matrix(c(TRUE, FALSE, TRUE,
             FALSE, TRUE, TRUE), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), c("C1", "C2", "C3"))),
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

  knn <- Matrix::Matrix(0, nrow = 3, ncol = 3, sparse = TRUE)
  rownames(knn) <- colnames(knn) <- c("C3", "C1", "C2")
  knn["C1", "C2"] <- 1
  knn["C3", "C2"] <- 1

  lcs <- IdentifyLogicConsensus(reo, knn_mat = knn, lr_db = lr_db, verbose = FALSE)
  expect_equal(unname(lcs["L_R"]), 0.5)

  bad_knn <- Matrix::Matrix(0, nrow = 2, ncol = 2, sparse = TRUE)
  rownames(bad_knn) <- colnames(bad_knn) <- c("C1", "C2")
  expect_error(
    IdentifyLogicConsensus(reo, knn_mat = bad_knn, lr_db = lr_db, verbose = FALSE),
    "contain all reo_mat column names"
  )
})

test_that("multimer logic requires all subunits", {
  logic_mat <- Matrix::Matrix(
    matrix(c(1, 1, 1,
             1, 0, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("A", "B"), paste0("C", 1:3))),
    sparse = TRUE
  )
  expect_equal(.resolve_complex_logic(c("A", "B"), logic_mat), c(TRUE, FALSE, TRUE))
  expect_equal(.resolve_complex_logic(c("A", "B", "C"), logic_mat), rep(FALSE, 3))
})

test_that("invalid arguments fail loudly", {
  lr_g <- all_lr_genes(lr_pairs_human)
  reo <- calc_REO_matrix(count_mat, lr_genes = lr_g, verbose = FALSE)
  scores <- score_lr_activity(reo, mode = "sender", aggregate = FALSE, verbose = FALSE)

  expect_error(score_lr_activity(reo, mode = "bad", verbose = FALSE))
  expect_error(rank_comm_cells(scores))

  lcs_list <- list(S1 = setNames(c(0.1, 0.0), c("A_B", "C_D")),
                   S2 = setNames(c(0.0), c("A_B")))
  expect_error(CompareLogicGroups(lcs_list, group_info = c(S1 = "Case"), verbose = FALSE))

  res <- CompareLogicGroups(
    list(S1 = setNames(c(0.1), "A_B"), S2 = setNames(c(0.0), "A_B")),
    group_info = c(S1 = "Case", S2 = "Ctrl"),
    verbose = FALSE
  )
  expect_error(filter_lcs(res, direction = "bad"))
})

test_that("REO anchors are computed from all genes by default, not only retained LR genes", {
  expr <- matrix(
    c(100, 100,
       90,  90,
       10,  10,
        9,   9),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("BG1", "BG2", "L", "R"), c("C1", "C2"))
  )
  reo <- calc_REO_matrix(expr, lr_genes = c("L", "R"),
                         rank_threshold = 0.5, verbose = FALSE)
  expect_equal(rownames(reo), c("L", "R"))
  expect_false(any(as.matrix(reo) != 0))
  expect_equal(attr(reo, "anchor_n_genes"), 4L)
})

test_that("positive_gate = any accepts any active positive modulator", {
  reo <- Matrix::Matrix(
    matrix(c(1, 1, 1,
             1, 1, 1,
             1, 0, 0,
             0, 1, 0), nrow = 4, byrow = TRUE,
           dimnames = list(c("L", "R", "AG1", "AG2"), paste0("C", 1:3))),
    sparse = TRUE
  )
  lr_db <- data.frame(
    lr_pair = "L_R", ligand = "L", receptor = "R",
    pathway = "test", annotation = "test", stringsAsFactors = FALSE
  )
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  lr_db$agonist_genes <- list(c("AG1", "AG2"))
  lr_db$antagonist_genes <- list(character(0))
  lr_db$co_A_receptor_genes <- list(character(0))
  lr_db$co_I_receptor_genes <- list(character(0))

  gate_any <- IdentifyLogicGateConsensus(
    reo, lr_db = lr_db, positive_gate = "any", negative_gate = "ignore",
    agonist_side = "receiver", verbose = FALSE
  )
  gate_all <- IdentifyLogicGateConsensus(
    reo, lr_db = lr_db, positive_gate = "all", negative_gate = "ignore",
    agonist_side = "receiver", verbose = FALSE
  )
  expect_equal(gate_any$base_lcs, 1)
  expect_equal(gate_any$regulated_lcs, 2 / 3)
  expect_equal(gate_all$regulated_lcs, 0)
})

test_that("CompareLogicGroups treats missing LCS as unavailable and preserves custom metadata", {
  lcs_list <- list(
    Case1 = setNames(0.2, "A_B"),
    Case2 = setNames(numeric(0), character(0)),
    Ctrl1 = setNames(0.0, "A_B"),
    Ctrl2 = setNames(0.0, "A_B")
  )
  groups <- c(Case1 = "Case", Case2 = "Case", Ctrl1 = "Ctrl", Ctrl2 = "Ctrl")
  lr_db <- data.frame(
    lr_pair = "A_B", ligand = "A", receptor = "B",
    pathway = "custom_pathway", annotation = "custom_db",
    stringsAsFactors = FALSE
  )
  res <- CompareLogicGroups(lcs_list, groups, lr_db = lr_db, verbose = FALSE)
  expect_equal(res$case_freq[res$lr_pair == "A_B"], 1)
  expect_equal(res$ctrl_freq[res$lr_pair == "A_B"], 0)
  expect_equal(res$n_case_avail[res$lr_pair == "A_B"], 1)
  expect_equal(res$pathway[res$lr_pair == "A_B"], "custom_pathway")
})

test_that("filter_lcs preserves valid data.frame structure and selected attributes", {
  res <- CompareLogicGroups(
    list(S1 = setNames(0.2, "A_B"), S2 = setNames(0.0, "A_B")),
    group_info = c(S1 = "Case", S2 = "Ctrl"), verbose = FALSE
  )
  attr(res, "custom_unused") <- "do not copy blindly"
  out <- filter_lcs(res, min_asymmetry = 0, max_fdr = 1, direction = "both")
  expect_true(inherits(out, "LogicCommResult"))
  expect_equal(nrow(out), length(attr(out, "row.names")))
  expect_equal(attr(out, "case_label"), attr(res, "case_label"))
  expect_equal(attr(out, "lcs_mat"), attr(res, "lcs_mat"))
  expect_null(attr(out, "custom_unused"))
})

test_that("KNN self loops are removed by default", {
  reo <- Matrix::Matrix(
    matrix(c(TRUE, FALSE,
             TRUE, FALSE), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), c("C1", "C2"))),
    sparse = TRUE
  )
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  knn <- Matrix::Matrix(0, nrow = 2, ncol = 2, sparse = TRUE)
  rownames(knn) <- colnames(knn) <- c("C1", "C2")
  knn["C1", "C1"] <- 1
  knn["C1", "C2"] <- 1

  lcs_drop <- IdentifyLogicConsensus(reo, knn_mat = knn, lr_db = lr_db,
                                     remove_self_edges = TRUE, verbose = FALSE)
  lcs_keep <- IdentifyLogicConsensus(reo, knn_mat = knn, lr_db = lr_db,
                                     remove_self_edges = FALSE, verbose = FALSE)
  expect_equal(unname(lcs_drop["L_R"]), 0)
  expect_equal(unname(lcs_keep["L_R"]), 0.5)
})

test_that("cell-type communication summaries recover directed A to B LCS and roles", {
  reo <- Matrix::Matrix(
    matrix(c(1, 1, 0, 0,
             0, 0, 1, 1,
             0, 1, 0, 0,
             0, 0, 1, 0), nrow = 4, byrow = TRUE,
           dimnames = list(c("L", "R", "X", "Y"), paste0("C", 1:4))),
    sparse = TRUE
  )
  labels <- setNames(c("A", "A", "B", "B"), paste0("C", 1:4))
  lr_db <- data.frame(
    lr_pair = c("L_R", "X_Y"),
    ligand = c("L", "X"),
    receptor = c("R", "Y"),
    pathway = c("p1", "p2"),
    annotation = "test",
    stringsAsFactors = FALSE
  )
  lr_db$ligand_genes <- list("L", "X")
  lr_db$receptor_genes <- list("R", "Y")

  knn <- Matrix::Matrix(0, nrow = 4, ncol = 4, sparse = TRUE)
  rownames(knn) <- colnames(knn) <- paste0("C", 1:4)
  knn["C1", "C3"] <- 1
  knn["C2", "C4"] <- 1
  knn["C3", "C1"] <- 1

  ct <- summarize_celltype_communication(
    reo, cell_labels = labels, knn_mat = knn, lr_db = lr_db,
    lcs_threshold = 0.5, min_edges = 1, verbose = FALSE
  )

  expect_true(inherits(ct, "LogicCommCellTypeComm"))
  row_ar <- ct$lr_table[ct$lr_table$sender_type == "A" &
                          ct$lr_table$receiver_type == "B" &
                          ct$lr_table$lr_pair == "L_R", ]
  expect_equal(row_ar$n_edges, 2)
  expect_equal(row_ar$n_active_edges, 2)
  expect_equal(row_ar$lcs, 1)

  ps <- ct$pair_summary[ct$pair_summary$sender_type == "A" &
                          ct$pair_summary$receiver_type == "B", ]
  expect_equal(ps$n_active_lr, 1)
  expect_equal(ps$sum_lcs, 1)
  expect_true(all(c("outdegree_score", "indegree_score", "betweenness_score", "information_score", "dominant_role") %in% names(ct$role_summary)))

  vec <- celltype_comm_to_lcs(ct, level = "celltype_lr", metric = "lcs")
  expect_true("A|B|L_R" %in% names(vec))
  expect_equal(unname(vec["A|B|L_R"]), 1)
})

test_that("hybrid cell-type scoring separates local, global, and denominator support", {
  reo <- Matrix::Matrix(
    matrix(c(1, 1, 0, 0,
             0, 0, 1, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), paste0("C", 1:4))),
    sparse = TRUE
  )
  labels <- setNames(c("A", "A", "B", "B"), paste0("C", 1:4))
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      pathway = "p1", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")

  knn <- Matrix::Matrix(0, nrow = 4, ncol = 4, sparse = TRUE)
  rownames(knn) <- colnames(knn) <- paste0("C", 1:4)
  knn["C1", "C2"] <- 1
  knn["C3", "C4"] <- 1

  ct_distal <- summarize_celltype_communication(
    reo, cell_labels = labels, knn_mat = knn, lr_db = lr_db,
    lcs_threshold = 0.5, min_edges = 1, min_active_edges = 1,
    verbose = FALSE
  )
  row_ab <- ct_distal$lr_table[ct_distal$lr_table$sender_type == "A" &
                                ct_distal$lr_table$receiver_type == "B" &
                                ct_distal$lr_table$lr_pair == "L_R", ]
  expect_equal(nrow(row_ab), 1)
  expect_true(row_ab$distal_candidate)
  expect_true(row_ab$active)
  expect_equal(row_ab$communication_range, "distal/endocrine")
  ps_ab <- ct_distal$pair_summary[ct_distal$pair_summary$sender_type == "A" &
                                    ct_distal$pair_summary$receiver_type == "B", ]
  expect_equal(ps_ab$n_distal, 1)
  expect_equal(ps_ab$n_local_active, 0)
  expect_equal(ps_ab$local_support_fraction_active, 0)
  expect_equal(ps_ab$communication_support_label, "global_only_candidate")

  ct_global <- summarize_celltype_communication(
    reo, cell_labels = labels, lr_db = lr_db, mode = "global",
    lcs_threshold = 0.5, min_edges = 1, min_active_edges = 1,
    verbose = FALSE
  )
  row_global <- ct_global$lr_table[ct_global$lr_table$sender_type == "A" &
                                    ct_global$lr_table$receiver_type == "B", ]
  expect_true(row_global$global_candidate_active)
  expect_true(row_global$active)
  expect_equal(ct_global$pair_summary$n_active_lr[ct_global$pair_summary$sender_type == "A" &
                                                    ct_global$pair_summary$receiver_type == "B"], 1)

  ct_min_edges <- summarize_celltype_communication(
    reo, cell_labels = labels, lr_db = lr_db, mode = "global",
    lcs_threshold = 0.5, min_edges = 5, min_active_edges = 1,
    verbose = FALSE
  )
  row_min_edges <- ct_min_edges$lr_table[ct_min_edges$lr_table$sender_type == "A" &
                                          ct_min_edges$lr_table$receiver_type == "B", ]
  expect_false(row_min_edges$global_candidate_active)
  expect_false(row_min_edges$active)
  expect_true(all(c("celltype_sizes", "adjacency_strength", "adjacency_count",
                    "adjacency_active_edge_support") %in% names(ct_min_edges)))
  expect_true(all(c("lcs_neighborhood", "lcs_global", "local_active",
                    "global_candidate_active", "distal_candidate", "candidate_active") %in% names(ct_min_edges$lr_table)))
})

test_that("cell-type plotting helpers return ggplot objects", {
  reo <- Matrix::Matrix(
    matrix(c(1, 1, 0, 0,
             0, 0, 1, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), paste0("C", 1:4))),
    sparse = TRUE
  )
  labels <- setNames(c("A", "A", "B", "B"), paste0("C", 1:4))
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      pathway = "p1", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  ct <- summarize_celltype_communication(
    reo, cell_labels = labels, lr_db = lr_db,
    mode = "global", lcs_threshold = 0.1, min_edges = 1, verbose = FALSE
  )
  expect_true(all(c("plot_communication_range_summary", "plot_lr_bubble_advanced") %in% getNamespaceExports("LogicComm")))
  expect_s3_class(plot_celltype_heatmap(ct), "ggplot")
  expect_s3_class(plot_celltype_roles(ct), "ggplot")
  expect_s3_class(plot_celltype_network(ct), "ggplot")
  expect_s3_class(plot_communication_range_summary(ct), "ggplot")
  expect_s3_class(plot_lr_bubble_by_celltype(ct, active_only = FALSE), "ggplot")
  expect_s3_class(plot_lr_bubble_advanced(ct), "ggplot")
  expect_s3_class(plot_pathway_heatmap(ct), "ggplot")
  expect_s3_class(plot_pathway_heatmap(ct, cluster_rows = TRUE, cluster_cols = TRUE), "ggplot")

  reo_loop <- Matrix::Matrix(
    matrix(c(1, 1,
             1, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), c("C1", "C2"))),
    sparse = TRUE
  )
  labels_loop <- setNames(c("A", "A"), c("C1", "C2"))
  ct_loop <- summarize_celltype_communication(
    reo_loop, cell_labels = labels_loop, lr_db = lr_db,
    mode = "global", lcs_threshold = 0.1, min_edges = 1, verbose = FALSE
  )
  p_loop <- plot_celltype_network(ct_loop, min_weight = 0, color_edges_by = "range")
  expect_s3_class(p_loop, "ggplot")
  expect_s3_class(ggplot2::ggplotGrob(p_loop), "gtable")
  expect_match(p_loop$labels$caption, "global-only")

  ex <- explain_celltype_interaction(ct, sender = "A", receiver = "B", lr_pair = "L_R")
  expect_true(is.list(ex))
  expect_true("evidence" %in% names(ex))
})

test_that("weighted graph LCS uses edge weights", {
  reo <- Matrix::Matrix(
    matrix(c(1, 1, 0,
             0, 1, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), c("C1", "C2", "C3"))),
    sparse = TRUE
  )
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  knn <- Matrix::Matrix(0, nrow = 3, ncol = 3, sparse = TRUE)
  rownames(knn) <- colnames(knn) <- c("C1", "C2", "C3")
  knn["C1", "C2"] <- 9
  knn["C3", "C2"] <- 1
  lcs_bin <- IdentifyLogicConsensus(reo, knn_mat = knn, lr_db = lr_db,
                                    edge_weight_mode = "binary", verbose = FALSE)
  lcs_w <- IdentifyLogicConsensus(reo, knn_mat = knn, lr_db = lr_db,
                                  edge_weight_mode = "weighted", verbose = FALSE)
  expect_equal(unname(lcs_bin["L_R"]), 0.5)
  expect_equal(unname(lcs_w["L_R"]), 0.9)
})

test_that("cell-type role summary exposes disambiguated count fields and low communication", {
  reo <- Matrix::Matrix(
    matrix(c(1, 0, 0, 0,
             0, 0, 1, 0), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), paste0("C", 1:4))),
    sparse = TRUE
  )
  labels <- setNames(c("A", "A", "B", "C"), paste0("C", 1:4))
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      pathway = "p", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  ct <- summarize_celltype_communication(
    reo, cell_labels = labels, lr_db = lr_db, mode = "global",
    lcs_threshold = 0.1, min_edges = 1, min_role_event_count = 2,
    verbose = FALSE
  )
  expect_true(all(c("outgoing_lr_event_count", "incoming_lr_event_count",
                    "outgoing_target_type_count", "incoming_source_type_count",
                    "outgoing_unique_lr_count", "role_confidence_label") %in% names(ct$role_summary)))
  expect_true(any(ct$role_summary$dominant_role == "Low-communication"))
})

test_that("cell-type role helpers are quiet for degenerate role scores", {
  adj <- matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE,
                dimnames = list(c("A", "B"), c("A", "B")))
  expect_no_warning(role_summary <- .communication_role_summary(adj, adj, min_role_event_count = 0))
  expect_true(is.data.frame(role_summary))
  expect_equal(.rescale01(c(NA_real_, NA_real_)), c(0, 0))
})

test_that("new cell-type visualization and uncertainty helpers return expected objects", {
  reo <- Matrix::Matrix(
    matrix(c(1, 1, 0, 0,
             0, 0, 1, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), paste0("C", 1:4))),
    sparse = TRUE
  )
  labels <- setNames(c("A", "A", "B", "B"), paste0("C", 1:4))
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      pathway = "p1", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  ct <- summarize_celltype_communication(
    reo, cell_labels = labels, lr_db = lr_db, mode = "global",
    lcs_threshold = 0.1, min_edges = 1, verbose = FALSE
  )
  expect_s3_class(plot_celltype_role_heatmap(ct), "ggplot")
  expect_s3_class(plot_celltype_role_radar(ct), "ggplot")
  expect_s3_class(plot_lr_activity_balance(ct, active_only = FALSE), "ggplot")
  expect_s3_class(plot_lr_evidence(ct, sender = "A", receiver = "B", lr_pair = "L_R"), "ggplot")
  boot <- bootstrap_celltype_communication(ct, n_boot = 5, seed = 1)
  expect_true(is.data.frame(boot))
  null_pair <- permute_celltype_communication(ct, reo_mat = reo, n_perm = 5, metric = "sum_lcs", seed = 1, verbose = FALSE)
  expect_true(is.data.frame(null_pair))
  expect_true(all(c("empirical_p", "n_null_nonmissing", "metric") %in% names(null_pair)))
  expect_equal(ct$params$edge_weight_mode, "binary")
  expect_warning(
    permute_celltype_communication(
      ct, reo_mat = reo, n_perm = 2, metric = "sum_lcs", seed = 1,
      edge_weight_mode = "weighted", verbose = FALSE
    ),
    "differs from ct_comm"
  )
  adaptive_messages <- character(0)
  adaptive_pair <- withCallingHandlers(
    permute_celltype_communication(
      ct, reo_mat = reo, n_perm = 2, adaptive = TRUE, adaptive_top_n = 1,
      adaptive_n_perm = 4, metric = "sum_lcs", seed = 1, verbose = TRUE
    ),
    message = function(m) {
      adaptive_messages <<- c(adaptive_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(is.data.frame(adaptive_pair))
  expect_true(any(grepl("Starting adaptive refinement", adaptive_messages)))

  reo_unlabeled <- Matrix::Matrix(
    matrix(c(1, 1, 0,
             0, 1, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), c("C1", "C2", "C3"))),
    sparse = TRUE
  )
  labels_unlabeled <- setNames(c("A", "B", NA_character_), c("C1", "C2", "C3"))
  ct_unlabeled <- summarize_celltype_communication(
    reo_unlabeled, cell_labels = labels_unlabeled, lr_db = lr_db, mode = "global",
    lcs_threshold = 0.1, min_edges = 1, verbose = FALSE
  )
  expect_equal(names(ct_unlabeled$cell_labels), c("C1", "C2"))
  null_unlabeled <- permute_celltype_communication(
    ct_comm = ct_unlabeled, reo_mat = reo_unlabeled, n_perm = 3,
    metric = "sum_lcs", seed = 1, verbose = FALSE
  )
  expect_true(is.data.frame(null_unlabeled))
  expect_equal(unique(null_unlabeled$metric), "sum_lcs")

  response_db <- data.frame(lr_pair = "L_R", response_genes = "R", stringsAsFactors = FALSE)
  resp <- score_receiver_response(ct, reo, response_db)
  expect_true("response_integrated_score" %in% names(resp))
})

test_that("publication extensions support spatial graphs, dynamics, and reports", {
  reo <- Matrix::Matrix(
    matrix(c(1, 1, 0, 0,
             0, 0, 1, 1), nrow = 2, byrow = TRUE,
           dimnames = list(c("L", "R"), paste0("C", 1:4))),
    sparse = TRUE
  )
  labels <- setNames(c("A", "A", "B", "B"), paste0("C", 1:4))
  coords <- data.frame(x = c(0, 0, 1, 1), y = c(0, 1, 0, 1), row.names = paste0("C", 1:4))
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      pathway = "p1", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")

  g <- build_spatial_graph(coords, k = 1, verbose = FALSE)
  expect_true(inherits(g, "sparseMatrix"))
  expect_equal(dim(g), c(4, 4))

  ct <- summarize_spatial_communication(
    reo, coords = coords, cell_labels = labels, lr_db = lr_db,
    k = 1, lcs_threshold = 0.1, min_edges = 1, verbose = FALSE
  )
  expect_true(inherits(ct, "LogicCommCellTypeComm"))
  expect_true("celltype_sizes" %in% names(ct))
  expect_s3_class(plot_spatial_logic(reo, coords, "L_R", lr_db = lr_db), "ggplot")
  expect_s3_class(plot_communication_qc(ct), "ggplot")
  expect_s3_class(plot_role_confidence(ct), "ggplot")

  findings <- summarize_communication_findings(ct, top_n = 2)
  expect_true(inherits(findings, "LogicCommFindings"))
  expect_true("distal_candidate_pairs" %in% names(findings))
  expect_false(any(findings$top_celltype_pairs$communication_support_label == "global_only_candidate"))
  tf <- tempfile(fileext = ".md")
  expect_true(file.exists(write_communication_report(ct, file = tf)))
  report_text <- paste(readLines(tf, warn = FALSE), collapse = "\n")
  expect_match(report_text, "Global-only distal candidate pairs")

  dyn <- summarize_communication_dynamics(
    reo, pseudotime = setNames(c(0.1, 0.2, 0.8, 0.9), paste0("C", 1:4)),
    cell_labels = labels, lr_db = lr_db, n_bins = 2,
    min_cells_per_bin = 1, mode = "global", lcs_threshold = 0.1,
    min_edges = 1, verbose = FALSE
  )
  expect_true(inherits(dyn, "LogicCommDynamics"))
  expect_s3_class(plot_communication_dynamics(dyn, level = "pair"), "ggplot")
})

test_that("sample-level GLM fits celltype communication features", {
  make_ct <- function(case_boost = FALSE) {
    reo <- Matrix::Matrix(
      matrix(c(1, 1, if (case_boost) 1 else 0, if (case_boost) 1 else 0,
               0, 0, 1, 1), nrow = 2, byrow = TRUE,
             dimnames = list(c("L", "R"), paste0("C", 1:4))),
      sparse = TRUE
    )
    labels <- setNames(c("A", "A", "B", "B"), paste0("C", 1:4))
    lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                        pathway = "p1", stringsAsFactors = FALSE)
    lr_db$ligand_genes <- list("L")
    lr_db$receptor_genes <- list("R")
    summarize_celltype_communication(
      reo, cell_labels = labels, lr_db = lr_db, mode = "global",
      lcs_threshold = 0.1, min_edges = 1, verbose = FALSE
    )
  }
  sample_ct <- list(Ctrl1 = make_ct(FALSE), Ctrl2 = make_ct(FALSE),
                    Case1 = make_ct(TRUE), Case2 = make_ct(TRUE))
  meta <- data.frame(group = factor(c("Ctrl", "Ctrl", "Case", "Case"), levels = c("Ctrl", "Case")),
                     row.names = names(sample_ct))
  fit <- fit_celltype_comm_glm(sample_ct, meta, design = ~ group, coef = "groupCase",
                               min_samples = 2, min_total_trials = 1, verbose = FALSE)
  expect_true(is.data.frame(fit))
  expect_true(all(c("feature", "estimate", "p_value", "fdr") %in% names(fit)))
  if (nrow(fit) > 0) expect_s3_class(plot_celltype_glm_volcano(fit), "ggplot")
})

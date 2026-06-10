test_that("logic validation adapters attach downstream evidence", {
  lr <- data.frame(
    sender_type = "A",
    receiver_type = "B",
    lr_pair = "L_R",
    ligand = "L",
    receptor = "R",
    pathway = "test",
    lcs = 0.5,
    active = TRUE,
    n_edges = 4,
    n_active_edges = 2,
    stringsAsFactors = FALSE
  )
  ct <- structure(
    list(
      lr_table = lr,
      cell_labels = c(C1 = "A", C2 = "B", C3 = "B")
    ),
    class = "LogicCommCellTypeComm"
  )
  reo <- Matrix::Matrix(
    matrix(c(1, 0, 1,
             0, 1, 1,
             0, 1, 0), nrow = 3, byrow = TRUE,
           dimnames = list(c("L", "R", "T1"), c("C1", "C2", "C3"))),
    sparse = TRUE
  )

  response_db <- data.frame(lr_pair = "L_R", response_genes = "T1", stringsAsFactors = FALSE)
  response <- logic_validate_receiver_response(ct, reo, response_db, attach = FALSE)
  expect_equal(response$validation_response_gene_count, 1L)
  expect_equal(response$validation_response_active_frac, 0.5)
  expect_equal(response$validation_response_score, 0.5)
  expect_equal(response$validation_response_integrated_lcs, 0.25)

  target_evidence <- data.frame(
    lr_pair = c("L_R", "L_R"),
    target_gene = c("T1", "T2"),
    target_score = c(0.8, 0.2),
    target_active_frac = c(1, 0.5),
    evidence_source = c("NicheNet", "NicheNet"),
    stringsAsFactors = FALSE
  )
  ct_target <- logic_validate_targets(ct, target_evidence)
  expect_equal(ct_target$lr_table$validation_target_n, 2L)
  expect_equal(ct_target$lr_table$validation_target_score_max, 0.8)
  expect_equal(ct_target$lr_table$validation_target_top_genes, "T1;T2")
  expect_equal(ct_target$validation$ligand_targets$validation_target_source, "NicheNet")

  tf_activity <- data.frame(
    receiver_type = c("B", "B"),
    tf = c("TF1", "TF2"),
    tf_activity_score = c(0.7, 0.3),
    tf_activity_delta = c(1.2, 0.4),
    tf_activity_fdr = c(0.01, 0.2),
    evidence_source = c("SCENIC", "SCENIC"),
    stringsAsFactors = FALSE
  )
  lr_tf_map <- data.frame(lr_pair = c("L_R", "L_R"), tf = c("TF1", "TF2"), stringsAsFactors = FALSE)
  ct_tf <- logic_validate_tf_switch(ct, tf_activity, lr_tf_map)
  expect_equal(ct_tf$lr_table$validation_tf_n, 2L)
  expect_equal(ct_tf$lr_table$validation_tf_score_max, 0.7)
  expect_equal(ct_tf$lr_table$validation_tf_switch_fdr_min, 0.01)
  expect_equal(ct_tf$lr_table$validation_tf_top, "TF1;TF2")

  external <- data.frame(
    sender_type = "A",
    receiver_type = "B",
    lr_pair = "L_R",
    evidence_score = 0.9,
    evidence_label = "protein_support",
    evidence_source = "proteomics",
    stringsAsFactors = FALSE
  )
  ct_ext <- logic_add_external_evidence(ct, external, evidence_type = "protein")
  expect_equal(ct_ext$lr_table$validation_external_n, 1L)
  expect_equal(ct_ext$lr_table$validation_external_score_mean, 0.9)
  expect_equal(ct_ext$lr_table$validation_external_type, "protein")
})

test_that("logic_grade_evidence grades active cell-type calls and validation evidence", {
  lr <- data.frame(
    sender_type = c("A", "A", "B"),
    receiver_type = c("B", "C", "A"),
    lr_pair = c("L_R", "L_X", "Y_R"),
    lcs = c(0.5, 0.4, 0),
    active = c(TRUE, TRUE, FALSE),
    validation_target_n = c(2L, 0L, 0L),
    n_edges = c(4, 4, 1),
    n_active_edges = c(2, 0, 0),
    stringsAsFactors = FALSE
  )
  ct <- structure(list(lr_table = lr), class = "LogicCommCellTypeComm")

  graded <- logic_grade_evidence(ct, min_edges = 2, min_active_edges = 1)
  # Row 1: active + validation -> A; row 2: active, no validation -> B;
  # row 3: not active -> insufficient (grade NA, tier 1).
  expect_equal(graded$lr_table$evidence_grade, c("A", "B", NA))
  expect_equal(graded$lr_table$evidence_tier, c(3L, 2L, 1L))
  expect_equal(graded$lr_table$evidence_components[1], "cell_type_coexpression;validation")
  expect_equal(graded$lr_table$evidence_components[2], "cell_type_coexpression")
  expect_equal(graded$lr_table$evidence_components[3], "insufficient_evidence")
  expect_match(graded$lr_table$interpretation_caution[3], "Edge opportunity count")
})

test_that("logic validation adapters reject malformed inputs", {
  ct <- structure(
    list(lr_table = data.frame(sender_type = "A", receiver_type = "B", lr_pair = "L_R", stringsAsFactors = FALSE)),
    class = "LogicCommCellTypeComm"
  )

  expect_error(
    logic_validate_targets(ct, data.frame(lr_pair = "L_R", target_gene = "T1")),
    "target_evidence is missing required columns"
  )
  expect_error(
    logic_validate_tf_switch(ct, data.frame(receiver_type = "B", tf = "TF1"), data.frame(lr_pair = "L_R")),
    "tf_activity is missing required columns"
  )
  expect_error(
    logic_add_external_evidence(ct, data.frame(unrelated = "x")),
    "evidence_table must share at least one identity column"
  )
})

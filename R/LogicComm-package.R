#' LogicComm: Logic-Based Cell-Cell Communication Analysis
#'
#' LogicComm identifies cell-cell communication by converting single-cell gene
#' expression into Relative Expression Ordering (REO) binary logic states, then
#' detecting ligand-receptor Logic Consensus Scores (LCS) across cells or cell
#' neighborhoods.
#'
#' Main workflow:
#' \enumerate{
#'   \item \code{\link{calc_REO_matrix}}: convert expression to a binary logic matrix.
#'   \item \code{\link{IdentifyLogicConsensus}}: compute per-sample LCS values.
#'   \item \code{\link{CompareLogicGroups}}: compare LCS between sample groups.
#'   \item \code{\link{summarize_celltype_communication}}: summarize cell-type resolved communication, receiver-response evidence, uncertainty diagnostics, and signaling roles.
#'   \item \code{\link{summarize_communication_dynamics}}: analyze pseudotime/time-bin communication.
#'   \item \code{\link{fit_celltype_comm_glm}}: fit sample-level differential communication models with covariates.
#'   \item \code{\link{run_multisample}}: run the full multi-sample pipeline.
#' }
#'
#' Visualization helpers include \code{\link{plot_lcs_bubble}},
#' \code{\link{plot_lcs_heatmap}}, \code{\link{plot_umap_logic}},
#' \code{\link{plot_celltype_heatmap}}, \code{\link{plot_celltype_roles}},
#' \code{\link{plot_celltype_role_heatmap}}, \code{\link{plot_celltype_network}},
#' \code{\link{plot_lr_bubble_by_celltype}}, \code{\link{plot_lr_evidence}},
#' \code{\link{plot_pathway_dominance}}, \code{\link{plot_celltype_role_dotplot}},
#' \code{\link{plot_specificity_stability}}, \code{\link{plot_celltype_network_publication}},
#' \code{\link{plot_communication_qc}}, and \code{\link{plot_celltype_glm_volcano}}.
#'
#' @keywords internal
"_PACKAGE"

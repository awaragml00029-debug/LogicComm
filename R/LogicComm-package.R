#' LogicComm: Logic-Based Cell-Cell Communication Analysis
#'
#' LogicComm identifies cell-cell communication by converting single-cell gene
#' expression into Relative Expression Ordering (REO) binary logic states, then
#' detecting ligand-receptor Logic Consensus Scores (LCS) at the cell-type
#' co-expression level.
#'
#' Main pipeline API:
#' \enumerate{
#'   \item \code{\link{logic_import_lrdb}} / \code{\link{logic_check_lrdb}}:
#'     bring a ligand-receptor database into LogicComm format.
#'   \item \code{\link{logic_prepare}}: convert expression to REO logic input.
#'   \item \code{\link{logic_score_lr}}: score ligand-receptor logic for one sample.
#'   \item \code{\link{logic_summarize_celltypes}}: summarize directed
#'     sender-cell-type to receiver-cell-type communication.
#'   \item \code{\link{logic_compare_groups}} or \code{\link{logic_run}}:
#'     compare sample groups or run the multi-sample workflow end to end.
#' }
#'
#' Advanced compatibility functions such as \code{\link{calc_REO_matrix}},
#' \code{\link{IdentifyLogicConsensus}}, and \code{\link{CompareLogicGroups}}
#' remain available for existing scripts, but new analyses should start from the
#' \code{logic_*} pipeline above.
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

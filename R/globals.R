# Avoid R CMD check notes for ggplot2 aesthetics and data-masked columns.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "lr_pair", "pathway", "log2fc_lcs", "asymmetry", "case_freq", "ctrl_freq",
    "UMAP1", "UMAP2", "status", "gene", "shift_score", "minus_log10_fdr",
    "lig_shift", "rec_shift", "x", "y", "xend", "yend", "value",
    "receiver_type", "sender_type", "hub_score", "sender_receiver_balance",
    "cell_type", "outdegree_score", "indegree_score", "dominant_role",
    "lr_pair_f", "pathway_plot", "color_val", "pathway", "celltype_pair",
    "role", "ligand_active_frac_sender", "receptor_active_frac_receiver",
    "component", "xval", "neglog10", "color_value", "label",
    ".data", "communication_range", "dominant_communication_range",
    "edge_color", "edge_linetype", "lcs", "axis_label", "Pathway",
    "Count", "Range", "metric", "feature", "to_label", "active_fraction",
    "pair_specificity", "total_lcs", "mean_lcs", "lr_pairs_human",
    "estimate", "top_pathway", "xm", "ym", "edge_text", "raw_score",
    "evidence", "sender_role_score", "receiver_role_score", "bin_order",
    "n_edges", "sum_lcs", "n_active_lr", "n_local_active", "n_distal_candidate",
    "mean_edge_support_fraction_active", "communication_support_label", "edge_support_alpha",
    "point_alpha", "size_val", "fraction", "neg_log10_fdr", "color_var", "role_confidence",
    "specificity_class", "participation_type", "share", "grp_key", "direction"
  ))
}

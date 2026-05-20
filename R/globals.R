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
    "component", "xval", "neglog10", "color_value", "label"
  ))
}

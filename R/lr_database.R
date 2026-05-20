#' Built-in human ligand-receptor pair database
#'
#' A curated data frame of 109 human ligand-receptor pairs. The seed set was
#' compiled from public ligand-receptor resources including CellChatDB,
#' CellPhoneDB, and NicheNet, then normalized into the compact LogicComm schema
#' used by the scoring functions.
#'
#' CellChatDB is a manually curated database of literature-supported ligand-
#' receptor interactions. In CellChat, the database is stored as a list with
#' components such as `interaction`, `geneInfo`, `complex`, and `cofactor`; the
#' human database is distributed as `CellChatDB.human.rda` in the CellChat
#' repository. LogicComm does not vendor the full CellChatDB object. Instead,
#' `lr_pairs_human` keeps a small, package-ready subset with explicit ligand and
#' receptor subunit list columns for REO/LCS scoring.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{lr_pair}{Interaction name, e.g. "VEGFA_KDR"}
#'   \item{ligand}{Ligand gene symbol (may be composite: "GeneA+GeneB")}
#'   \item{receptor}{Receptor gene symbol (may be composite: "GeneC+GeneD")}
#'   \item{ligand_genes}{List column of individual ligand subunit gene symbols}
#'   \item{receptor_genes}{List column of individual receptor subunit gene symbols}
#'   \item{pathway}{Signaling pathway name}
#'   \item{annotation}{Functional category}
#' }
#'
#' @examples
#' data(lr_pairs_human)
#' head(lr_pairs_human[, c("lr_pair","ligand","receptor","pathway")])
#' all_lr_genes(lr_pairs_human)[1:10]
"lr_pairs_human"

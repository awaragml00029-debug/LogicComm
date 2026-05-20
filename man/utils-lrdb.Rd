\name{all_lr_genes}
\alias{all_lr_genes}
\alias{as_logiccomm_lr_db_from_cellchat}
\title{Ligand-receptor database utilities}
\description{Utilities for extracting genes from LR databases and converting CellChatDB-like objects.}
\usage{all_lr_genes(lr_db, include_modulators = FALSE)
as_logiccomm_lr_db_from_cellchat(cellchat_db, pathway_col = "pathway_name",
  annotation_col = "annotation", interaction_name_col = "interaction_name",
  source_label = "CellChatDB", make_unique = TRUE,
  include_modulators = TRUE)}
\arguments{
\item{lr_db}{LogicComm LR database.}
\item{include_modulators}{Whether to include modulator list columns when present.}
\item{cellchat_db}{CellChatDB-like list with an \code{interaction} data frame and optional \code{complex} table.}
\item{pathway_col}{Interaction column used as pathway.}
\item{annotation_col}{Interaction column used as annotation.}
\item{interaction_name_col}{Interaction column used as pair identifier.}
\item{source_label}{Annotation label used when \code{annotation_col} is absent.}
\item{make_unique}{Whether to make duplicate LR pair names unique.}
}
\value{\code{all_lr_genes()} returns a character vector; \code{as_logiccomm_lr_db_from_cellchat()} returns a data frame.}
\examples{data(lr_pairs_human)
all_lr_genes(lr_pairs_human)[1:5]
}

\name{logic_validation}
\alias{logic_validation}
\alias{logic_validate_receiver_response}
\alias{logic_validate_targets}
\alias{logic_validate_tf_switch}
\alias{logic_add_external_evidence}
\alias{logic_grade_evidence}
\title{LogicComm validation adapter API}
\description{
Attach downstream or external evidence to LogicComm cell-type communication events.
These helpers keep communication inference separate from validation while making it
straightforward to connect receiver-response genes, ligand-target evidence,
TF/regulon activity, protein/spatial/perturbation evidence, and evidence grades.
}
\usage{
logic_validate_receiver_response(ct_comm, reo_mat, response_db,
  response_mode = c("any", "all"), attach = TRUE)

logic_validate_targets(ct_comm, target_evidence,
  by = c("lr_pair", "ligand"), attach = TRUE)

logic_validate_tf_switch(ct_comm, tf_activity, lr_tf_map, attach = TRUE)

logic_add_external_evidence(ct_comm, evidence_table, evidence_type = NULL,
  attach = TRUE)

logic_grade_evidence(ct_comm, use_validation = TRUE, min_edges = NULL,
  min_active_edges = NULL)
}
\arguments{
\item{ct_comm}{A \code{LogicCommCellTypeComm} object, usually returned by
\code{logic_summarize_celltypes()} or \code{summarize_celltype_communication()}.}
\item{reo_mat}{Binary REO matrix or \code{LogicCommREOResult}.}
\item{response_db}{Data frame with \code{lr_pair} and \code{response_genes}.
Response genes may be separated by \code{;}, \code{,}, \code{|}, or \code{+}.}
\item{response_mode}{How to combine response genes: \code{"any"} marks a
receiver cell active when any response gene is active; \code{"all"} requires all
response genes.}
\item{attach}{If \code{TRUE}, return the updated communication object. If
\code{FALSE}, return only the validation table.}
\item{target_evidence}{NicheNet-style ligand-target evidence table with
\code{target_gene}, \code{target_score}, and either \code{lr_pair} or
\code{ligand}. Optional \code{sender_type}, \code{receiver_type},
\code{target_active_frac}, and \code{evidence_source} columns add specificity and
summary fields.}
\item{by}{Whether target evidence is matched by \code{"lr_pair"} or
\code{"ligand"}.}
\item{tf_activity}{Receiver TF or regulon activity table with columns
\code{receiver_type}, \code{tf}, and \code{tf_activity_score}. Optional columns
include \code{tf_activity_delta}, \code{tf_activity_fdr}, and
\code{evidence_source}.}
\item{lr_tf_map}{Mapping table connecting communication events to TFs. Must
contain \code{tf} and share at least one of \code{lr_pair}, \code{pathway}, or
\code{ligand} with \code{ct_comm$lr_table}.}
\item{evidence_table}{External evidence table sharing at least one identity
column with \code{ct_comm$lr_table}, such as \code{sender_type},
\code{receiver_type}, \code{lr_pair}, \code{ligand}, \code{receptor}, or
\code{pathway}. Optional columns include \code{evidence_score},
\code{evidence_label}, \code{evidence_source}, and \code{evidence_type}.}
\item{evidence_type}{Optional evidence type to assign when
\code{evidence_table} does not already contain \code{evidence_type}.}
\item{use_validation}{Whether validation columns should contribute to evidence
grades.}
\item{min_edges}{Optional edge-opportunity threshold used to add caution text.}
\item{min_active_edges}{Optional active-edge threshold used to add caution text.}
}
\details{
The validation API is intentionally adapter-oriented. LogicComm first reports
communication evidence from REO/LCS and graph or global cell-type logic. These
functions then attach independent support: receiver read-receipt genes, predicted
or measured ligand targets, TF/regulon switches, and external evidence such as
protein, spatial, perturbation, or literature support.

\code{logic_grade_evidence()} adds publication-facing interpretation columns.
Grade \code{A} indicates local LogicComm evidence plus validation support;
\code{B} indicates local LogicComm evidence; \code{C} indicates global or distal
candidate potential; \code{D} indicates weak or limited active evidence. Grades
are evidence summaries, not causal proof.
}
\value{
\code{logic_validate_receiver_response()}, \code{logic_validate_targets()},
\code{logic_validate_tf_switch()}, and \code{logic_add_external_evidence()}
return either an updated \code{LogicCommCellTypeComm} object or a validation data
frame when \code{attach = FALSE}. \code{logic_grade_evidence()} returns an
updated \code{LogicCommCellTypeComm} object with evidence-grade columns in
\code{lr_table}.
}
\examples{
\dontrun{
ct_comm <- logic_summarize_celltypes(reo, cell_labels = celltypes, lr_db = lr_db)

response_db <- data.frame(
  lr_pair = "TGFB1_TGFBR1",
  response_genes = "SMAD7;SERPINE1",
  stringsAsFactors = FALSE
)
ct_comm <- logic_validate_receiver_response(ct_comm, reo, response_db)

nichenet_targets <- data.frame(
  lr_pair = "TGFB1_TGFBR1",
  target_gene = c("SMAD7", "SERPINE1"),
  target_score = c(0.8, 0.7),
  stringsAsFactors = FALSE
)
ct_comm <- logic_validate_targets(ct_comm, nichenet_targets)

ct_comm <- logic_grade_evidence(ct_comm)
}
}

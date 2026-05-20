\name{LogicComm-package}
\alias{LogicComm-package}
\alias{LogicComm}
\docType{package}
\title{Logic-Based Cell-Cell Communication Analysis}
\description{LogicComm converts single-cell expression into Relative Expression Ordering (REO) binary logic states and scores ligand-receptor communication with Logic Consensus Scores (LCS), including weighted graph scoring, cell-type communication, receiver-response evidence, uncertainty diagnostics, spatial/pseudotime extensions, sample-level GLM modeling, publication reports, and signaling-role summaries.}
\details{The main workflow is \code{calc_REO_matrix()}, \code{IdentifyLogicConsensus()}, and \code{CompareLogicGroups()} or the one-stop \code{run_multisample()} wrapper. Use \code{summarize_celltype_communication()} for celltypeA -> celltypeB interaction counts, strengths, pathway summaries, and sender/receiver/mediator/influencer role analysis. Use \code{bootstrap_celltype_communication()}, \code{permute_celltype_communication()}, and \code{sensitivity_REO_threshold()} for robustness checks. Use \code{build_spatial_graph()} and \code{summarize_spatial_communication()} for spatial neighborhoods, \code{summarize_communication_dynamics()} for pseudotime/time bins, \code{fit_celltype_comm_glm()} for sample-level differential modeling, and \code{write_communication_report()} for markdown reporting.}
\keyword{package}

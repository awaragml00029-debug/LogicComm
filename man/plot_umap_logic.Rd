\name{plot_umap_logic}
\alias{plot_umap_logic}
\title{UMAP plot colored by LR logic activity}
\description{Color cells by sender/receiver/both/inactive state for a selected LR pair.}
\usage{plot_umap_logic(reo_mat, seurat_obj, lr_pair, lr_db = lr_pairs_human,
  pt_size = 0.8, colors = NULL)}
\arguments{
\item{reo_mat}{Binary REO matrix.}
\item{seurat_obj}{Seurat object with UMAP coordinates.}
\item{lr_pair}{LR pair name.}
\item{lr_db}{Ligand-receptor database.}
\item{pt_size}{Point size.}
\item{colors}{Optional named colors.}
}
\value{A ggplot object.}

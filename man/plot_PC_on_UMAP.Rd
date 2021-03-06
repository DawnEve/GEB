% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCA.R
\name{plot_PC_on_UMAP}
\alias{plot_PC_on_UMAP}
\title{Plot cell PC score as color on UMAP plot}
\usage{
plot_PC_on_UMAP(sce, nPCs = 16, ncol = 3)
}
\arguments{
\item{sce}{Seurat object}

\item{nPCs}{an integer, total number of PC to plot;}
}
\value{
a ggplot2 merge obj, by cowplot::plot_grid
}
\description{
Plot cell PC score as color on UMAP plot
}
\examples{
plot_PC_on_UMAP(pbmc, ncol=4)
plot_PC_on_UMAP(pbmc, 9)
}

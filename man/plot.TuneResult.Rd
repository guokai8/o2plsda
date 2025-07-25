% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.TuneResult}
\alias{plot.TuneResult}
\title{Enhanced plot method for TuneResult objects}
\usage{
\method{plot}{TuneResult}(x, type = "heatmap", metric = "mean_score", ...)
}
\arguments{
\item{x}{A TuneResult object}

\item{type}{Character. Type of plot: "heatmap", "line", "comparison", "optimal"}

\item{metric}{Character. Performance metric to visualize}

\item{...}{Additional plotting parameters}
}
\description{
Enhanced plot method for TuneResult objects
}
\examples{
# Example with tuning results
\dontrun{
# Assuming you have a TuneResult object from tune_o2pls()
tune_results <- tune_o2pls(X, Y, nc_range = 1:3, nx_range = 0:2, ny_range = 0:2)

# Heatmap of results
plot(tune_results, type = "heatmap")

# Line plot showing trends
plot(tune_results, type = "line")

# Comparison across parameters
plot(tune_results, type = "comparison")
}

}

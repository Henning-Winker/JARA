% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jara_utils.R
\name{mixed.trend}
\alias{mixed.trend}
\title{mixed.trend()}
\usage{
mixed.trend(
  jara,
  run = "joint",
  refyr = FALSE,
  type = c("mu", "pr")[1],
  threshold.id = 1
)
}
\arguments{
\item{jara}{fit_jara output}

\item{run}{name qualifier of data.frame}

\item{refyr}{sets index relative to a reference year}

\item{type}{c("mu","pr"), note the probability is only available of "mixed.trends"}

\item{threshold.id}{option to select threshold location from list if multiple were specified}
}
\value{
data.frame
}
\description{
Function to compile combined index
}
\author{
Henning Winker (JRC-EC)
}

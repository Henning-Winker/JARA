% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_jara.R
\name{build_jara}
\alias{build_jara}
\title{build_jara function}
\usage{
build_jara(
  I = NULL,
  se = NULL,
  assessment = "Unnamed",
  scenario = "s1",
  model.type = c("relative", "census")[1],
  mixed.scale = c("geomean", "mean", "median")[1],
  GL = NULL,
  start.year = NA,
  end.year = NA,
  fixed.obsE = NULL,
  variance.weighting = "equal",
  sigma.proc.fixed = FALSE,
  proc.pen = NULL,
  timeblock = FALSE,
  proj.mod = c("theta10", "logistic", "theta.value")[1],
  pk.prior = c(0.25, 0.1),
  pk.yr = NULL,
  pk.i = NULL,
  proj.r = c("all", "GL1", "year")[1],
  proj.yrs.user = NULL,
  proj.stoch = FALSE,
  backcast.na = NULL,
  silent = FALSE
)
}
\arguments{
\item{I}{abundance indices/counts,  require data.frame(year, I.1,I.2,...,I.N)}

\item{se}{optional log standard error (CV) time series,requires data.frame(year, se.1,se.2,...,se.N)}

\item{assessment}{= "species.X",}

\item{scenario}{= "s1",}

\item{model.type}{abundance data type c("relative","census","mixed.trends")}

\item{mixed.scale}{choice if "mixed.trends" are presented by the geomean, mean or median}

\item{GL}{Generation length. If not specified by the user (recommended) then the default is GL = floor((n.years-3)/3)}

\item{start.year}{subsetting option for start year}

\item{end.year}{subsetting option for end year}

\item{fixed.obsE}{minimum plausible observation error (fixed)}

\item{variance.weighting}{option "equal", "model", "fixed" or user assinged vector e.g. c(1,1,2)}

\item{sigma.proc.fixed}{option to fix the process error e.g = 0.1 (default FALSE)}

\item{proc.pen}{advanced user setting to penalize extreme process error deviations
New break point settings}

\item{timeblock}{option to set a time block for mean.r a specified year (default FALSE)}

\item{proj.mod}{specified by theta c("theta10","logistic")}

\item{pk.prior}{lognormal prior of population to K c(mu,lod.sd,yr) envoked during projections}

\item{pk.yr}{reference year for the depletion prior. NULL ratio to max pred. abundance across all years}

\item{pk.i}{option to specify specify different depletion levels for colonies in the census model}

\item{proj.r}{# rate of change for c("all","GL1", years), "all" is default}

\item{proj.yrs.user}{option to overwrite GL and costomize projections for forecasting}

\item{proj.stoch}{allows for projections with process error c(TRUE, FALSE), FALSE is default}

\item{backcast.na}{option to remove NAs bootstrap when backcasting by specifying NA}

\item{silent}{option to put notifications on siltent}
}
\value{
List to be used as data input to JARA JAGS model
}
\description{
Creates a data list with JABBA input and settings to be passed to fit_jara()
}
\examples{
data(jaradata)
inp = jrdat$Afr_penguin
buildjr <- build_jara(I = inp$I,GL=9,model.type="census")
jrplot_indices(buildjr)
}
\author{
Henning Winker
}

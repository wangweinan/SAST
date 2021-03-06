% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAST.R
\name{SAST}
\alias{SAST}
\title{Perform SAST Online FDR procedure}
\usage{
SAST(z_score, alpha = 0.05, conservative = TRUE, init = 200, h)
}
\arguments{
\item{z_score}{standardized observations}

\item{alpha}{nominal FDR level, default is 0.05}

\item{conservative}{logical value. Default is TRUE. When TRUE, run the conservative version of SAST by estimating non-null proportion as zero in constructing the test statistics. Otherwise, estimate varying proportion.}

\item{init}{initial burn in period of observations at the beginning of provided z scores. Default is 200.}

\item{h}{bandwidth choice for covariates}
}
\description{
This function implements the main SAST procedure in Gang et al. (2020). It takes a vector of z-scores (standardized observations)
and return the decision made at nominal FDR level alpha (default is 0.05).
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNM_function.R
\name{CN.R}
\alias{CN.R}
\title{CN Clustering}
\usage{
CN.R(data_x, it_num, inf, sup)
}
\arguments{
\item{data_x}{sample from GMM}

\item{it_num}{number of iterations}

\item{inf}{inf bound of the range for searching for parameters}

\item{sup}{sup bound of the range for searching for parameters}
}
\value{
the clustering parameters and proportions
}
\description{
Implement constrained Newton method to estimate the parameters and mixture proportions.
}
\examples{
\dontrun{
data_x = c(rnorm(500,0),rnorm(500,2));
M = CN.R(data_x, it_num = 20, inf = -1.5, sup = 3.5);M
}
}

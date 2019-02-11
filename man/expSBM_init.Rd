\name{expSBM_init}
\alias{expSBM_init}
\title{expSBM_init}
\description{
Initialization step for the variational expectation maximization algorithm. 
}
\usage{
expSBM_init(edgelist, K, soft = TRUE) 
}
\arguments{
  \item{edgelist}{A matrix with 4 columns: on the first column the sender node, on the second the receiver, on the third either a one or zero to indicate whether it is an interaction or a non-interaction respectively, on the fourth the corresponding exponential length.}
  \item{K}{Number of latent groups.}
  \item{soft}{\code{TRUE} or \code{FALSE} to indicate whether the function should return a soft clustering or not.}
}
\value{
A \code{NxK} matrix indicating the partitioning of the nodes.
}
\examples{
set.seed(1)
data(high_school)
K <- 4
lambda_init <- rep(1/K,K)
expSBM_init(high_school$edgelist, K)$Z
}



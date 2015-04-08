\name{cheby-package}
\alias{cheby-package}
\alias{cheby}
\docType{package}
\title{
Computes polynomial approximations of arbitrary functions
}
\description{
Package to implement Chebychev and shape-preserving Chebychev approximations in one and two dimensions.
}
\details{
\tabular{ll}{
Package: \tab cheby\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-24\cr
License: \tab GPL-3\cr
}

Generates polynomial approximations of arbitrary functions in one and two dimensions, whilst preserving slope, concavity and higher derivatives wherever specified.

\code{d1.poly} is the main routine for approximating one-dimensional funcitons.  \code{sp1.poly} does the same but preserves the sign of an arbitrary number of derivatives.  2-Dimensional Chebychev polynomials are also available with \code{dn.poly}.  Higher order Chebychev and shape-preserving approximations to follow later.

Comments and suggestions are gratefully received by the author.
}
\author{
Philip Barrett <pobarrett@uchicago.edu>
}
\references{
Judd, Kenneth L (1998) Numerical Methods in Economics
}
\keyword{ package }
\seealso{
\code{\link{d1.poly}}, \code{\link{sp1.poly}}, \code{\link{dn.poly}}
}
\examples{
## Compute basic approximations to natural logarithm
RR <- d1.poly( log, c(0,4), 6, 10 )
SS <- sp1.poly( log, c(0,4), 6, 10, n.shape=c(5,10), sign.deriv=c(1,-1), solver='NLOPT_LD_SLSQP' )
pp <-  seq( 0, 4, length.out=100 )
plot( pp, sapply(pp, RR), lwd=2, col=2, type='l' )
lines( pp, sapply(pp, log), lwd=2, col=1 )
lines( pp, sapply(pp, SS), lwd=2, col=4 )
legend( 'bottomright', c( 'log', 'Order 6 polynomial approx', 'Order 6 shape-preserving polynomial approx' ), lwd=2, col=c(1,2,4), bty='n' )
}

\name{sp1.poly}
\alias{sp1.poly}
\title{Shape-preserving polynomial approximation}
\usage{
sp1.poly(fn, range, iOrder, iPts, fn.opts = NULL, fn.vals = NULL,
  grid = NULL, n.shape = 0, sign.deriv = NULL, x0 = NULL,
  solver = "NLOPT_LD_SLSQP", tol = 1e-06, poly.return = FALSE,
  quiet = FALSE)
}
\arguments{
  \item{fn}{a function \eqn{f(x)} or \eqn{f(x,\beta)} for
  \eqn{\beta} a list of function parameters.  If the
  latter, must be coded with second argument a list names
  \code{opts}, i.e. \code{fn <- function( x, opts )}}

  \item{range}{the range of the approximation.}

  \item{iOrder}{the order of the polynomial approximation.}

  \item{iPts}{the number of points at which the
  approximation is computed. Must be at least as large as
  \code{iOrder}.}

  \item{fn.opts}{(optional) options passed to \code{fn}}

  \item{fn.vals}{the values of \code{fn} on \code{grid}.
  Useful if \code{fn} is very slow to evaluate.}

  \item{grid}{(optional) the grid on which the function is
  to be approximated.}

  \item{n.shape}{a vector of the number of shape-preserving
  points for each order of differentiation.  For example,
  to specify the slope at 5 Chebychev points and the
  concavity at 10, use \code{n.shape=c( 5, 10 )}}

  \item{sign.deriv}{a vector of signs {+1, 0, -1 } defining
  the sign of each derivative.  For example, for a concave
  approximation with positive slope sign.deriv=c(1,-1).}

  \item{x0}{initial guess of the polynomial approximation.
  Default is the standard Chebychev approximation.}

  \item{solver}{the \code{nlopt} solver to use in computing
  the best fit.  Default is \code{NLOPT_LD_SLSQP}.}

  \item{tol}{tolerance for solver convergence.  Default is
  \code{1e-06}}

  \item{poly.return}{If \code{TRUE}, returns extra details
  about the approximation.}
}
\value{
A function which approximates \code{fn}.  If
\code{poly.return=TRUE}, returns also the polynomial
description ovver [-1,1]
}
\description{
Approximates the function \code{fn} using shape-preserving
polynomials evaluated at the points in grid.
}
\examples{
base <- d1.poly( log, c(0,4), 6, 10, details=TRUE )
sp.compare <- sp1.poly(  log, c(0,4), 6, 10, poly.return=TRUE )
sp.flat.x0 <- sp1.poly(  log, c(0,4), 6, 10, x0=c(1,1,1,1,1,1,1), poly.return=TRUE )
print( base$poly - sp.compare$poly )
print( base$poly - sp.flat.x0$poly )
   # Comparison without using shape-preserving methods
sp.concave <- sp1.poly( log, c(0,4), 6, 10, n.shape=c(5,10), sign.deriv=c(1,-1) )
pp <-  seq( 0, 4, length.out=100 )
plot( pp, sapply(pp, base$fn), lwd=2, col=2, type='l' )
lines( pp, sapply(pp, log), lwd=2, col=1 )
lines( pp, sapply(pp, sp.concave), lwd=2, col=4 )
legend( 'bottomright', c( 'log', 'Order 6 polynomial approx',
     'Order 6 shape-preserving polynomial approx' ), lwd=2, col=c(1,2,4), bty='n' )
   # Compare the Chebychev and shape-preserving approximations
}
\references{
\code{nloptr} documentation
\link{http://cran.r-project.org/web/packages/nloptr/nloptr.pdf}
}
\seealso{
\code{\link{d1.poly}}
}


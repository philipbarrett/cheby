\name{d1.poly}
\alias{d1.poly}
\title{1-dimensional Chebychev approximation}
\usage{
d1.poly(fn, range, iOrder, iPts, fn.opts = NULL, fn.vals = NULL,
  grid = NULL, details = FALSE)
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

  \item{details}{If \code{TRUE}, returns extra details
  about the approximation.}
}
\value{
A function which approximates the input fn over the
interval \code{range}. If \code{details=TRUE}, also
includes the polynomial desciption over [-1,1], as well as
the approximation errors
}
\description{
Standard Chebychev approximation of an arbitrary function.
}
\examples{
cube <- function( x, opts ) opts$A * x^3
approx <- d1.poly( cube, c(-4,2), 4, 20, fn.opts=list(A=2) )
sapply( c(-3, -2, 0, .5 ), function( x ) abs( approx(x) - 2 * x ^ 3 ) )
}
\seealso{
\code{\link{sp1.poly}}
}


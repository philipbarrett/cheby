#####################################################################################
# d1.R - computes the basic Chebychev polynomial approximation of a one-dimensional #
# function                                                                          #
# Philip Barrett, Chicago                                                           #
# Created: 23dec2013                                                                #
#####################################################################################

#' Calculate grid of approximating points
#' 
#' Computes either a Uniform or Chebychev grid of collocation nodes
#' @param vRange the vector of the range over which the nodes are computed
#' @param iPts number of colocations points
#' @param Either \code{Chebychev} or \code{Uniform}
#' @export
d1.grid <- function( vRange, iPts, stMethod='Chebychev' ){
  if( stMethod=='Chebychev' ){
    zz <- sapply( 1:iPts, function(i) return( - cos( ( 2 * i - 1 ) / ( 2 * iPts ) * pi ) ) )
    return( ( zz + 1 ) * diff( vRange ) / 2 + vRange [ 1 ] )
  }       # Chebychev interpolation nodes
  if ( stMethod=='Uniform' ) 
    return( seq( from=vRange[1], to=vRange[2], length.out=iPts ) )
  # The uniform grid
}

#' @export
d1.normalize <- function( x, range )
  return( 2 * ( x - min( range ) ) / diff( range ) - 1 )
# Converts [a,b] to [-1,1]


#' 1-dimensional Chebychev approximation
#' 
#' Standard Chebychev approximation of an arbitrary function.
#' 
#' @param fn a function \eqn{f(x)} or \eqn{f(x,\beta)} for \eqn{\beta} a list of
#'   function parameters.  If the latter, must be coded with second argument a 
#'   list names \code{opts}, i.e. \code{fn <- function( x, opts )}
#' @param range the range of the approximation.
#' @param iOrder the order of the polynomial approximation.
#' @param iPts the number of points at which the approximation is computed. Must
#'   be at least as large as \code{iOrder}.
#' @param fn.opts (optional) options passed to \code{fn}
#' @param fn.vals the values of \code{fn} on \code{grid}.  Useful if \code{fn} 
#'   is very slow to evaluate.
#' @param grid (optional) the grid on which the function is to be approximated.
#' @param details If \code{TRUE}, returns extra details about the approximation.
#'   
#' @return A function which approximates the input fn over the interval 
#'   \code{range}. If \code{details=TRUE}, return is a list with entries 
#'   \code{fn, poly, fn.deriv, poly.deriv, residuals}, which are, respectively, 
#'   the approximating function, the polynomial desciption over [-1,1], the 
#'   derivative of the approximation, the polynomial desciption of the
#'   derivative, and the approximation errors.
#' @seealso \code{\link{sp1.poly}}
#' @examples
#' cube <- function( x, opts ) opts$A * x^3
#' approx <- d1.poly( cube, c(-4,2), 4, 20, fn.opts=list(A=2) )
#' sapply( c(-3, -2, 0, .5 ), function( x ) abs( approx(x) - 2 * x ^ 3 ) )
#' @export
d1.poly <- function( fn, range, iOrder, iPts, fn.opts=NULL, fn.vals=NULL, 
                     grid=NULL, details=FALSE ){
  
  stopifnot( iOrder <= iPts )
        # Check number of points is at least as big as the grid
  if ( is.null( grid ) ) grid <- d1.grid( range, iPts )
        # Computes the Chebychev grid if not supplied
  if( is.null( fn.vals ) ) fn.vals <- 
    if( is.null( fn.opts ) ) sapply( grid, fn) else sapply( grid, fn, opts=fn.opts )
        # Computes the function over the grid if not already submitted, passing
        # function options if required
  zz <- sapply( grid, d1.normalize, range=range )
        # Compress the grid into [-1,1]
  basis.poly <- chebyshev.t.polynomials( iOrder, normalized=FALSE )
        # A list of polynomial objects
  basis.vals <- as.matrix( sapply( zz, function(z.pt) 
                            sapply( basis.poly, function( pol ) as.function( pol )( z.pt ) ) ) )
        # The values of the basis on the grid
  basis.reg <- lm( fn.vals ~ 0 + t( basis.vals ) )
        # The regression of the function evaluation points on the basis grid
        # (zero eliminates intercept)
  poly <- do.call( 'sum', lapply( 1:(iOrder+1), function(i) basis.reg$coef[i] * basis.poly[[i]] ) )
          # Construct a polynomial  representation of the function
  fn.poly <- function( x ) as.function(poly)( d1.normalize( x, range ) ) 
        # The function returning the value of the polynomial over the range
  if( !details ) return( fn.poly )
  
  poly.deriv <- deriv( poly )
  f.deriv <- function( x ) 2 / ( diff( range ) ) * as.function(poly.deriv)( d1.normalize( x, range ) ) 
  return( list( fn=fn.poly, poly=poly, fn.deriv=f.deriv, deriv.poly=poly.deriv, residuals=basis.reg$residuals ) )
  
}

d1.chebyconvert <- function( vCheby ){
# Converts a vector of Chebychev polynomial coefficients into an ordinary polynomial
  
  iOrder <- length( vCheby ) - 1
        # Need to correct for zero-order term
  basis.poly <- chebyshev.t.polynomials( iOrder, normalized=FALSE )
        # Generate the basis polynomials
  poly <- do.call( 'sum', lapply( 1:(iOrder+1), function(i) vCheby[i] * basis.poly[[i]] ) )
        # Convert
  return(poly)
}

d1.fn <- function( x, range, poly ){
  # Computes the value of poly(x) when the range of x is given by range, rather
  # than by [-1,1]
  return( as.function(poly)( d1.normalize( x, range ) ) )
}
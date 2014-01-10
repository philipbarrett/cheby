#####################################################################################
# dn.R - computes the basic Chebychev polynomial approximation of a multi-          #
# dimensional function                                                              #
# Philip Barrett, Chicago                                                           #
# Created: 08jan2013                                                                #
#####################################################################################

#' Mutli-dimensional Chebychev approximation
#' 
#' Standard Chebychev approximation of an arbitrary function. [INCOMPLETE -
#' CURRENTLY AIMING ONLY FOR 2-DIMENSIONAL FUNCTIONS]
#' 
#' @param fn a function \eqn{f(x_1,..., x_n)} or \eqn{f(x_1, ..., x_n,\beta)} 
#'   for \eqn{\beta} a list of function parameters.  If the latter, must be 
#'   coded with second argument a list names \code{opts}, i.e. \code{fn <- 
#'   function( x, opts )}
#' @param range the range of the approximation, given as a list of vectors. Eg 
#'   for a 3-dimensional function approximated over [1,2] * [-1,2] * [0,4], 
#'   would be \code{range = list( c(1,2), c(-1,2), c(0,4) )}
#' @param iOrder the vector of orders of the polynomial approximation. Eg. to 
#'   approximate using polynomials of order 3 and 5 in the 1st and 2nd 
#'   dimensions respectively, would be \code{iOrder=c(5,6)}
#' @param iPts the vector of number of points at which the approximation is 
#'   computed. Must be at least as large as \code{iOrder} (element-by-element).
#' @param fn.opts (optional) options passed to \code{fn} [NOT YET FUNCTIONAL]
#' @param fn.vals the values of \code{fn} on \code{grid}.  Useful if \code{fn} 
#'   is very slow to evaluate.
#' @param grid (optional) the grid on which the function is to be approximated. 
#'   Should be submitted as a list of vectors for the grids in each dimension.
#' @param details If \code{TRUE}, returns extra details about the approximation.
#'   
#' @return A function which approximates the input fn over the box defined by 
#'   \code{range}. If \code{details=TRUE}, also includes the polynomial 
#'   desciption over [-1,1], as well as the approximation errors
#' @seealso \code{\link{d1.poly}}
#' @examples
#' # tbc
#' @export
dn.poly <- function( fn, range, iOrder, iPts, fn.opts=NULL, fn.vals=NULL, 
                     grid=NULL, details=FALSE ){
  
  iDim <- length( iOrder )
  stopifnot( length( range ) == iDim, all( iOrder <= iPts ) )
      # Error checking
  if ( is.null( grid ) ) grid <- mapply( d1.grid, range, iPts )
      # Computes the Chebychev grids in each dimension (if not already supplied)
  if( is.null( fn.vals ) ) fn.vals <- 
    if( is.null( fn.opts ) ) 
      do.call( multi.outer, c( fn, grid ) )
    else stop( 'Function options not yet available' ) # sapply( grid, fn, opts=fn.opts )
        # Computes the function over the grid if not already submitted, passing
        # function options if required
  fn.vec <- as.vector( t( fn.vals ) )
        # The vector of function values
  zz <- mapply( d1.grid, rep( list(c(-1,1)), iDim ), iPts )
        # Compress the grid into [-1,1]
  basis.poly <- lapply( iOrder, chebyshev.t.polynomials, normalized=FALSE )
        # A list of (list of) polynomial objects
  basis.vals <- lapply( 1:iDim, function(i) as.matrix( sapply( zz[[i]], function(z.pt) 
                      sapply( basis.poly[[i]], function( pol ) as.function( pol )( z.pt ) ) ) ) )
        # The values of the basis on the grid
  
  # FROM HERE ON 2 DIMENSIONS ONLY
  basis.vec <- t( do.call( 'cbind', lapply( 1:iPts[1], function(i) apply( basis.vals[[2]], MARGIN=2, 
                    function(col2) as.vector( outer( basis.vals[[1]][,i], col2 ) ) ) ) ) )
        # A matrix values where each row is the value of the cross-product of
        # the polynomials at each grid point
  basis.reg <- lm( fn.vec ~ 0 + basis.vec )
        # The regression of the function evaluation points on the basis grid
        # (zero eliminates intercept)
  
  
  # AFTER HERE MY BRAIN FAILED
  poly <- do.call( 'sum', lapply( 1:(iOrder+1), function(i) basis.reg$coef[i] * basis.poly[[i]] ) )
        # Construct a polynomial  representation of the function
  fn.poly <- function( x ) as.function(poly)( d1.normalize( x, range ) ) 
        # The function returning the value of the polynomial over the range
  if( !details ) return( fn.poly )
  return( list( fn=fn.poly, poly=poly, residuals=basis.reg$residuals ) )
  
}

# These two functions create a multi-dimensional version of outer
# Taken from stackoverflow solution
# http://stackoverflow.com/questions/6192848/how-to-generalize-outer-to-n-dimensions
list_args <- Vectorize( function(a,b) c( as.list(a), as.list(b) ), 
                        SIMPLIFY = FALSE)
make_args_mtx <- function( alist )
  Reduce(function(x, y) outer(x, y, list_args), alist)
multi.outer <- function(f, ... ) {
  args <- make_args_mtx(list(...))
  apply(args, 1:length(dim(args)), function(a) do.call(f, a[[1]] ) )
 }
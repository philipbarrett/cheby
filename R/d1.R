#####################################################################################
# d1.R - computes the basic Chebychev polynomial approximation of a one-dimensional #
# function                                                                          #
# Philip Barrett, Chicago                                                           #
# Created: 23dec2013                                                                #
#####################################################################################

d1.grid <- function( vRange, iPts, stMethod='Chebychev' ){
# Computes the grid of approximating points
  
  if( stMethod=='Chebychev' ){
    zz <- sapply( 1:iPts, function(i) return( - cos( ( 2 * i - 1 ) / ( 2 * iPts ) * pi ) ) )
    return( ( zz + 1 ) * diff( vRange ) / 2 + vRange [ 1 ] )
  }       # Chebychev interpolation nodes
  if ( stMethod=='Uniform' ) 
    return( seq( from=vRange[1], to=vRange[2], length.out=iPts ) )
  # The uniform grid
}

d1.normalize <- function( x, range )
  return( 2 * ( x - min( range ) ) / diff( range ) - 1 )
# Converts [a,b] to [-1,1]

d1.poly <- function( fn, range, iOrder, iPts, fn.opts=NULL, fn.vals=NULL, grid=NULL, details=FALSE ){
# Approximates the function fn using Chebychev polynomials evaluated at the
# points in grid, returns a chebychev polynomial object
  
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
  return( list( fn=fn.poly, poly=poly, residuals=basis.reg$residuals, 
                std.err=summary(basis.reg)$coefficients[, 2 ] ) )
  
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
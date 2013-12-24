####################################################################################
# uncertaintyShocks.R - R file containing function definitions for computing value #
# functions for NGM planner's problem when the variance of productivity is         #
# (potentially) stochastic                                                         #
# Philip Barrett, Chicago                                                          #
# Created: 06nov2013                                                               #
####################################################################################

sp1.ssr <- function( grid, poly, fn.vals, range ){
# Returns the mean square of residual of the polynomial approximation
  poly.vals <- sapply( grid, d1.fn, range=range, poly=poly )
        # The values of the polynomial evaluated on the grid
  return( sum( ( poly.vals - fn.vals ) ^ 2 ) )
}

sp1.ssr.grad <- function( grid, poly, fn.vals, range ){
# Returns the vector of derivatives of sp1.ssr
  poly.vals <- sapply( grid, d1.fn, range=range, poly=poly )
        # The values of the polynomial evaluated on the grid
  coeffs.deriv <- sapply( 0:(length(poly)-1), function(i) sapply( grid, d1.normalize, range=range) ^ i )
        # The derivative of the polynomial wrt each coefficient is the power of
        # the value at that point.  Eg, if y = b_0 + b_1 x + b_2 x^2, then the 
        # Derivative wrt b is ( 1, x_1, x_2 )
  return( c( 2 * ( poly.vals - fn.vals ) %*% coeffs.deriv ) )
}

sp1.deriv.order <- function( grid, poly, range, order, return.poly=FALSE ){
# Computes, at a grid of points, the order-th derivative of a Chebychev
# polynomial over a finite range
  poly.deriv <- list( poly )
        # List structure for differentiation
  for (i in 1:order) poly.deriv <- polynomial.derivatives( poly.deriv )
        # Repeatedly differentiate until we have the order-th order differential
  poly.deriv <- poly.deriv[[1]]
        # Urgh.  Need to unlist without killing structure
  if( return.poly ) return( poly.deriv )
        # Return the polynomial form of the derivative
  fn.deriv <- function( X ) 
    return( ( 2 / diff( range ) ) ^ order * d1.fn( X, range, poly.deriv ) )
        # The gradient function.  We need the diff(range) to correct for the
        # fact that the approximation is on range, not [-1,1]
  return( sapply( grid, fn.deriv ) )
}

sp1.deriv.order.grad.fd <- function( grid, poly, range, order ){
# Finite difference method to differentiate sp1.deriv.order w.r.t. the individual
# terms of the polynomial
  
  poly.order <- length(poly)
        # The order+1 of the polynomial
  base <- sp1.deriv.order( grid, poly, range, order )
        # The reference value
  increment.matrix <- diag( poly.order ) * 1e-06
        # The matrix of increments for the FD calculation
  fd.mat <- sapply( 1:poly.order, function(i) ( sp1.deriv.order( grid, poly+polynomial( increment.matrix[i, ] ), 
                                                                    range, order ) - base ) / 1e-06 )
        # The matrix of first-differences
  return( fd.mat )
#   return( as.vector( t( fd.mat ) ) )
}

sp1.poly <- function( fn, range, n.order, n.points, fn.opts=NULL, fn.vals=NULL, grid=NULL, 
                            n.shape=0, sign.deriv=NULL, x0=NULL, solver='NLOPT_LD_SLSQP', 
                            tol=1e-06, poly.return=FALSE ){
# Approximates the function fn using shape-preserving Chebychev polynomials 
# evaluated at the points in grid, returns a chebychev polynomial object.  Shape
# constraints are given by: a vector of numbers, representing the number of 
# shape constrants of each order, n.shape; plus a vector of signs { +1, 0, -1 } 
# defining the sign of each derivative, sign.deriv.  Eg. for a concave
# approximation with positive slope sign.deriv=c(1,-1), and if the slope and
# curvature are evlauted at 5 and 10 points respectively, n.shape=c( 5, 10 )
  
  # 0. Set up
  if ( is.null( x0 ) ) x0 <- as.vector( d1.poly( fn, range, n.order, n.points, 
                                                 fn.opts, fn.vals, details= TRUE )$poly )
  ub <- Inf * rep( 1, n.order + 1 )
  lb <- -Inf * rep( 1, n.order + 1 )
  if ( is.null( grid ) ) grid <- d1.grid( range, n.points )
        # Computes the Chebychev grid if not supplied
  if( is.null( fn.vals ) ) fn.vals <- 
    if( is.null( fn.opts) ) sapply( grid, fn ) else sapply( grid, fn, opts=opts ) 
        # Computes the function over the grid if not already submitted
  
  # 1. Define optimization functions
  eval_f <- function( x ) return( sp1.ssr( grid, polynomial( x ), fn.vals, range ) )
  eval_grad_f <- function( x ) return( sp1.ssr.grad( grid, polynomial( x ), fn.vals, range ) )
  if ( sum( n.shape) != 0 ){
    shape.grids <- lapply( n.shape, function(iII) d1.grid( range, iII ) )
    eval_g <- function( x ) return( unlist( lapply( 1:length(n.shape), 
                function(i) - sign.deriv[i] * sp1.deriv.order( shape.grids[[i]], polynomial(x), range, i ) ) ) )
    eval_jac_g <- function( x ) return( do.call( 'rbind', lapply( 1:length(n.shape), 
                    function(i) - sign.deriv[i] * sp1.deriv.order.grad.fd( shape.grids[[i]], 
                                                                        polynomial(x), range, i ) ) ) )
  }

  # 2. Solve the optimization problem
  
  #* 2.1 The unconstrained problem *#
  if ( sum( n.shape) == 0 ){
    nlopt.opts <- list('algorithm'=solver, "xtol_rel"=tol)
    optimize <- nloptr(x0 = x0, eval_f = eval_f, eval_grad_f = eval_grad_f, lb = lb,
                        ub = ub, opts = nlopt.opts )
  }
  #* 2.2 With shape-preserving constraints *#
  else{
      nlopt.opts <- list('algorithm'=solver, "xtol_rel"=tol)
      optimize <- nloptr(x0 = x0, eval_f = eval_f, eval_grad_f = eval_grad_f, lb = lb, ub = ub, 
                         eval_g_ineq = eval_g, eval_jac_g_ineq = eval_jac_g, 
                         opts = nlopt.opts )
  }
  
  #* 3. Check the status of the solver and return *#
  if( optimize$status > 0)
    message( 'Sucess: nlopt terminated after ', optimize$iterations, ' iterations \n',
             optimize$message )
  if( optimize$status < 0)
    stop('nlopt did not solve OK')
  
  fn.out <- function( x ) as.function(polynomial( optimize$solution ) )( d1.normalize( x, range ) )
        # Create the funciton from the polynomial
  if( poly.return)
    return( list( poly = polynomial( optimize$solution ), fn=fn.out ) )
  else
    return( fn.out )
}

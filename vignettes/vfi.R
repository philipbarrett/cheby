#####################################################################################
# vfi.R - Value function iteration for the neoclassical growth model                #
#                                                                                   #
# Philip Barrett, Chicago                                                           #
# Created: 26apr2015                                                                #
#####################################################################################

library(nloptr)

util <- function( cons, sigma ){
# Computes period utility
  if( sigma == 1 ) return( log ( cons ) )
  return( ( cons ^ ( 1 - sigma ) ) / ( 1 - sigma ) )
}

optim.T <- function( kap, vf, alpha, betta, sigma, delta, k.range, pol=FALSE ){
# Computes the value of the maximization problem given a continuation function
# vf when capital is kap
  
  #* i. Specifics of the problem *#
  lb <- c( 0, k.range[1] )
  ub <- c( kap ^ alpha + ( 1 - delta ) * kap, k.range[2] )
      # Range of next-period capital choices
  x0 <- rep( .5 * ( kap ^ alpha + ( 1 - delta ) * kap ), 2 )
      # Initial guess that satsfies the constraint
  
  #* ii. Translate the functions above into the objective functions *#
  eval_f <- function( cons.kap )
    return( - ( util( cons.kap[1], sigma ) + betta * vf$fn( cons.kap[2] ) ) )
      # The maximand
  eval_grad_f <- function( cons.kap ) 
    return( c( - cons.kap[1] ^ (- sigma ), - betta * vf$fn.deriv( cons.kap[2] ) ) )
      # The gradient of the maximand
  eval_g_eq <- function( cons.kap )
    return( sum( cons.kap ) - ( kap ^ alpha + ( 1 - delta ) * kap ) )
      # The law of motion for capital
  eval_jac_g_eq <- function( cons.kap )
    return( c( 1, 1 ) )
      # The constraint jacobian
  
  #* iii. Solve the optimization problem *#
  nlopt.opts <- list('algorithm'='NLOPT_LD_SLSQP', "xtol_rel"=1.0e-8)
  optimize <- nloptr(x0 = x0, eval_f = eval_f, eval_grad_f = eval_grad_f, lb = lb,
                     ub = ub, eval_g_eq=eval_g_eq, eval_jac_g_eq=eval_jac_g_eq, 
                     opts = nlopt.opts )
  
  #* iv. Check the status of the solver and return *#
  if ( optimize$status < 0 )
    stop('NLOPT did not solve OK')
  if( pol==TRUE ) 
    return( optimize$solution )
  return( -optimize$objective )
}

vf.update <- function( vf, alpha, betta, sigma, delta, k.range, order, sp=FALSE ){
# Returns a new Chebychev approximation to the value function
  optim.T.params <- function( kap )
    return( optim.T( kap, vf, alpha, betta, sigma, delta, k.range ) )
      # The operator T with current parameters substituted in
  if( sp ){
    out <- sp1.poly( fn=optim.T.params, range=k.range, iOrder=order,
                     iPts=4*order, details=TRUE, n.shape=c(2, 2 * order ),
                     sign.deriv=c(1,-1), quiet=TRUE )
  }else{
    out <- d1.poly( fn=optim.T.params, range=k.range, iOrder=order, iPts=order+1, 
                    details=TRUE )
  }
      # The polynomial approximation to the operator
  return( out )
}

vf.iterate <- function( vf, alpha, betta, sigma, delta, order, 
                        max.diff, max.it, sp=FALSE ){
# Solves the neoclassical growth model via value function iteration
  
  k.range <- c( 1e-04, delta ^ ( 1 / ( alpha - 1 ) ) )
      # Range of k
  diff <- 2 * max.diff
  it <- 0
      # Initialize reporting variables
  test.grid <- seq( k.range[1], k.range[2], length.out=100 )
  vf.grid <- sapply( test.grid, vf$fn )
      # The grid of points to compute the change in value functions
  while( diff > max.diff & it < max.it ){
    vf.new <- vf.update( vf, alpha, betta, sigma, delta, k.range, order, sp )
    vf.grid.new <- sapply( test.grid, vf.new$fn )
    diff <- max( abs( vf.grid - vf.grid.new ) )
        # Compute the difference between the two value functions
    it <- it+1
    vf <- vf.new
    vf.grid <- vf.grid.new
        # assign new variables to old
  }
  optim.T.pols <- function( kap )
    return( optim.T( kap, vf, alpha, betta, sigma, delta, k.range, TRUE ) )
      # The operator T with current parameters substituted in.  Selects policies
  pols <- rbind( test.grid, sapply( test.grid, optim.T.pols ) )
      # The policy functions
  return( list( vf=vf, pols=pols, diff=diff, it=it, k.range=k.range ) )
}

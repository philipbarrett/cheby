#######################################################################################
# Test functions for d1 source file                                                   #
# Philip Barrett, Chicago                                                             #
# Created: 30jul2013                                                                  #
#######################################################################################


#### 1. ASSIGNMENT ####

test.grid <- function(){
# Testing d1.grid

  # Test that the central point of the Chebychev grid is always the midpoint of
  # the range
  grid <- d1.grid( c( 0,1), 3 )
  checkEquals( grid, .5 * c( cos( 5 * pi / 6 ) + 1, 1, cos( 1 * pi / 6 ) + 1 ) )
  grid <- d1.grid( c(4,8), 9 )
  checkEquals( grid[5], 6 )
  checkEquals( length( grid ), 9 )

  # Check that the uniform grid is as it should be
  grid <- d1.grid( c(4,8), 5, stMethod='Uniform' )
  checkEquals( grid, 4:8 )
}

test.normalize <- function(){
# Testing d1.normalize
  
  checkEquals( d1.normalize( 1, c( 0, 2 ) ), 0 )
  checkEquals( d1.normalize( 0, c( 0, 2 ) ), -1 )
  checkEquals( d1.normalize( 2, c( 0, 2 ) ), 1 )
  checkEquals( d1.normalize( 3, c( 0, 4 ) ), .5 )
}
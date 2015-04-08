#######################################################################################
# Test suite for cheby package                                                        #
# Philip Barrett, Chicago                                                             #
# Created: 24dec2013                                                                  #
#######################################################################################

## Now run the unit tests if RUnit is available
if(require("RUnit", quietly=TRUE)) {
  
  # 1. Run the test suite
#   require("RUnit")
  cat("\nRunning unit tests\n")
        # Screen updating
  test.suite <- defineTestSuite( 'cheby', dirs = file.path('inst/unitTests') )
        # Define the suite
  test.results <- runTestSuite( test.suite )
        # Record the results
  
  # 2. Report the results
  printTextProtocol( test.results, showDetails=FALSE )
        # Print the results  
  printTextProtocol( test.results, showDetails=TRUE,
                     fileName="tests/report.txt" )
  
  # 3. Stop if there are errors
  tmp <- getErrors( test.results )
  if(tmp$nFail > 0 | tmp$nErr > 0) {
    stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
               ", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
  }
  
} else {
  warning("cannot run unit tests -- package RUnit is not available")
}
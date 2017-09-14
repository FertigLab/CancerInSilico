#define CATCH_CONFIG_RUNNER
#include "Tests/catch.h"

// [[Rcpp::export]]
int run_catch_unit_tests()
{
    Catch::Session session;
    int numFailed = session.run();
    return (numFailed < 0xFF ? numFailed : 0xFF);
}
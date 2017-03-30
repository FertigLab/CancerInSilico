// [[Rcpp::depends(BH)]]

#include <boost/random.hpp>
#include <stdint.h>

#include "Random.h"

static boost::random::mt19937 rng;

void Random::setSeed(uint32_t newSeed)
{
    rng.seed(newSeed);
}

double Random::normal(double mean, double var)
{
    boost::random::normal_distribution<double> dist(mean, var);
    return dist(rng);
}

int Random::uniformInt(int a, int b)
{
    boost::random::uniform_int_distribution<> dist(a,b);
    return dist(rng);
}

double Random::uniform(double a, double b)
{
    boost::random::uniform_distribution<> dist(a,b);
    return dist(rng);
}


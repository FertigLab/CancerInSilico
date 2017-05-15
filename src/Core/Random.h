#ifndef CIS_RANDOM_H
#define CIS_RANDOM_H

// wrapper class for random number generator

#include <stdint.h>

namespace Random
{
    void setSeed(uint32_t);

    int uniformInt(int, int);
    double uniform(double, double);
    double normal(double, double);
}

#endif

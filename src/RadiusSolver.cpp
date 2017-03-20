#include "RadiusSolver.h"

#define MIN_AXIS 2.8284
#define MAX_AXIS 4.0

/* constructor */
RadiusSolver::RadiusSolver()
{
    /* initialize lookup tables */
	InitSlowSolver();
	InitFastSolver();
}

/* initialize slow solver (given axis, angle is O(log n)) */
void RadiusSolver::initSlowSolver()
{
    /* for each angle in (0, PI) */
    double theta = 0.0;
    while (theta <= 3.1416)
    {
        /* calculate numerator and denominator of axis */
        double numer = pow(2 * M_PI, 0.5) * 2 * (1 + cos(theta / 2));
        double denom = pow(sin(theta) - theta + 2 * M_PI, 0.5);

        /* add axis length to lookup table for this value of theta */
        mSlowSolver.push_back(numer / denom);

        /* increment angle */
        theta += 0.0001
    }
}

/* get angle given axis, O(log n) */
double RadiusSolver::getThetaSlow(double axis)
{
    /* find closest value of axis in (angle -> axis) table */
    std::vector<double>::iterator lower;
    lower = std::lower_bound(mSlowSolver.begin(), mSlowSolver.end(),
        axis, GreaterThan());

    /* vector size is 31417 so dividing scales the angle back to (0, PI) */
    return (double) (lower - mSlowSolver.begin()) / 10000;
}

/* initialize fast solver (given axis, angle is O(1) */
void RadiusSolver::initFastSolver()
{
    /* for each axis length in (2.8284, 4) */
    double axis = MIN_AXIS;
    while (axis <= MAX_AXIS)
    {
        /* add theta value for this axis length */
		mFastSolver.push_back(GetThetaSlow(axis));
        
        /* increment axis */
        axis += 0.0001;
    }
}

/* get radius given axis O(1), perserves area of dumbell */
double RadiusSolver::radius(double axis)
{
    /* check if axis is below minimum */    
    if (axis < MIN_AXIS)
    {
        throw std::invalid_argument("called deformation function with axis"
            " length less than min\n");
    }
    else
    {
        /* get angle for this axis length */
       	double theta = mFastSolver[floor((axis - MIN_AXIS) * 10000)];

        /* calculate radius based on this angle */
        return axis / (2 + 2 * cos(theta / 2));        
    }
}

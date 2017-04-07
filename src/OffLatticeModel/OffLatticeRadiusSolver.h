#ifndef CIS_OFF_LATTICE_RADIUS_SOLVER_H
#define CIS_OFF_LATTICE_RADIUS_SOLVER_H

#include <vector>

// used for sorting/finding
struct GreaterThan
{
    bool operator()(double a, double b) const
    {
        return a > b;
    }
};

/* This class manages the calculation of radius values after the axis
   has been deformed. A dumbell shaped cell has 3 variables that are
   dependent on each other. Axis, Angle, Radius always have the same
   relationship no matter the stage of division the cell is in. This
   class provides the new Radius after the Axis has been changed. The
   calculation in done by calculating the value of the Angle given the
   Axis, and the Radius given the Angle */
class OffLatticeRadiusSolver
{
private:

    /* lookup tables for radius-axis values */
	std::vector<double> mSlowSolver;
	std::vector<double> mFastSolver;

    /* initialize slow solver table (given axis, angle is O(log n)) */
    void initSlowSolver();

    /* initialize fast solver table (given axis, angle is O(1)) */
    void initFastSolver();

    /* get angle given axis, O(log n) */
    double getThetaSlow(double axis);

public:

    /* constructor */
    OffLatticeRadiusSolver();

    /* get radius given axis in O(1), perserves area of dumbell */
    double radius(double) const;
};

#endif

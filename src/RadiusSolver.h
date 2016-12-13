#ifndef RADIUS_SOLVER_H
#define RADIUS_SOLVER_H

#include <vector>

struct GreaterThan {

    bool operator()(double a, double b) const {
        
        return a > b;

    }

};

class RadiusSolver {

private:

    /* lookup tables for radius-axis values */
	std::vector<double> mSlowSolver;
	std::vector<double> mFastSolver;

    /* initialize lookup tables for radius-axis values */
	void InitSlowSolver();
	void InitFastSolver();

    /* hash axis length for fast lookup table */
	int HashAxisLength(double);

public:

    /* constructor */
    RadiusSolver();

    /* return radius given axis length, perserves area of dumbell */
    double GetRadius(double);

};

#endif

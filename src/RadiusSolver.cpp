#include "RadiusSolver.h"

RadiusSolver::RadiusSolver()
{
	InitSlowSolver();
	InitFastSolver();
}

void RadiusSolver::InitSlowSolver()
{
    double theta, numer, denom;

    for (int i = 0; i <= 31400; ++i)
    {
        theta = (double) i / 10000;
        numer = pow(2 * M_PI, 0.5) * 2 * (1 + cos(theta / 2));
        denom = pow(sin(theta) - theta + 2 * M_PI, 0.5);
        mSlowSolver.push_back(numer / denom);
    }
}

void RadiusSolver::InitFastSolver() {

	double theta;

	for (int i = 0; i <= 11716; ++i) {

		theta = GetThetaSlow(2 * pow(2,0.5) + (double) i / 10000);
		mFastSolver.push_back(theta);

	}

}

double RadiusSolver::GetThetaSlow(double axis) {

    std::vector<double>::iterator lower;
    lower = std::lower_bound(mSlowSolver.begin(), mSlowSolver.end(),
                                axis, GreaterThan());

    return (double) (lower - mSlowSolver.begin()) / 10000;

}

double RadiusSolver::GetRadius(double axis) {

    if (axis < 2 * pow(2,0.5)) {

        throw std::invalid_argument("called deformation function with axis"
            " length less than min\n");

    } else {

       	double theta = mFastSolver[HashAxisLength(axis)];
        return axis / (2 + 2 * cos(theta / 2));        

    }

}

int RadiusSolver::HashAxisLength(double axis) {

	return floor((axis - 2 * pow(2,0.5)) * 10000);

}


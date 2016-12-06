#include <cmath>
#include <Rcpp.h>
#include <algorithm>
#include <exception>

#include "Parameters.h"

Parameters::Parameters(double maxRad, Rcpp::List Rparams) {

    mParams = Rparams;
    mMaxRadius = maxRad;
	InitializeRadiusSolver();

    StoreTimeIncrement();
    StoreUpdateParameters();
    StoreGrowthDistribution();

}

void Parameters::StoreTimeIncrement() {

    Rcpp::NumericVector cycleDist = mParams["cycleLengthDist"];
    double minCycle = Rcpp::min(cycleDist);

    double t1 = delta() / (4 * nG() * (4 - pow(2, 0.5)));
    double t2 = delta() * (minCycle - 1) / (8 * nG() * (pow(2, 0.5) - 1));

    mParams["timeIncrement"] = std::min(t1, t2);

}

void Parameters::StoreUpdateParameters() {

    mParams["maxDeform"] = 2 * timeIncrement() * nG() * (4 - pow(2, 0.5));
    mParams["maxTranslation"] = delta() / 2;
    mParams["maxRotate"] = acos((16 + pow(delta(), 2) - 4 * delta()) / 16);

}

void Parameters::StoreGrowthDistribution() {

    Rcpp::NumericVector cycleDist = mParams["cycleLengthDist"];

    for (unsigned int i = 0; i < cycleDist.size(); ++i) {

        double numer = 2 * (pow(2, 0.5) - 1) * timeIncrement() * nG();
        mGrowthDist.push_back(numer / cycleDist[i]);

    }

}

void Parameters::InitializeRadiusSolver() {

	InitSlowSolver();
	InitFastSolver();

}

//hardcoded for max radius = sqrt(2)
void Parameters::InitSlowSolver() {

    double theta, numer, denom;

    for (int i = 0; i <= 31400; ++i) {

        theta = (double) i / 10000;
        numer = pow(2 * M_PI, 0.5) * 2 * (1 + cos(theta / 2));
        denom = pow(sin(theta) - theta + 2 * M_PI, 0.5);
        mSlowSolver.push_back(numer / denom);

    }

}

//m_slow_solver sorted from high to low
double Parameters::GetThetaSlow(double axis) {

    std::vector<double>::iterator lower;
    lower = std::lower_bound(mSlowSolver.begin(), mSlowSolver.end(),
                                axis, GreaterThan());

    return (double) (lower - mSlowSolver.begin()) / 10000;

}

void Parameters::InitFastSolver() {

	double theta;

	for (int i = 0; i <= 11716; ++i) {

		theta = GetThetaSlow(2 * pow(2,0.5) + (double) i / 10000);
		mFastSolver.push_back(theta);

	}

}

int Parameters::HashAxisLength(double axis) {

	return floor((axis - 2 * pow(2,0.5)) * 10000);

}

double Parameters::GetRadius(double axis) {

    if (axis < 2 * pow(2,0.5)) {

        throw std::invalid_argument("called deformation function with axis length less than max\n");

    } else {

       	double theta = mFastSolver[HashAxisLength(axis)];
        return axis / (2 + 2 * cos(theta / 2));        

    }

}

double Parameters::GetRandomGrowthRate() {

	return mGrowthDist[floor(R::runif(0, mGrowthDist.size()))];

}

double Parameters::GetDrugEffect(double growthRate) {

    Rcpp::Function de = mParams["drugEffect"];
    double cycleLength = 1 + 2 * (pow(2, 0.5) - 1) * timeIncrement() * nG() / growthRate;
    return Rcpp::as<double>(de(cycleLength));

}



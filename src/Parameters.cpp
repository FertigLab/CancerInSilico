#include <cmath>
#include <Rcpp.h>
#include <algorithm>

#include "ExceptionHandling.h"
#include "Parameters.h"

#define MIN_AXIS_LEN 1.1716

void Parameters::InitializeRadiusSolver() {

	InitSlowSolver();
	InitFastSolver();

}

//hardcoded for max radius = sqrt(2)
void Parameters::InitSlowSolver() {

    double theta, num, denom;

    for (int i = 0; i <= 31400; ++i) {

        theta = (double) i / 10000;
        num = 64 * M_PI * pow((1 + cos(theta / 2)), 2);
        denom = sin(theta) - theta + 2 * M_PI;
        m_slow_solver.push_back(pow(num / denom, 0.5));

    }

}

void Parameters::InitFastSolver() {

	double theta;

	for (int i = 0; i <= 11716; ++i) {

		theta = GetThetaSlow(2 * pow(2,0.5) + (double) i / 10000);
		m_fast_solver.push_back(theta);

	}

}

double Parameters::GetThetaSlow(double axis_len) {

    int imin = 0, imax = 31400;
    int imid;

    while (imin + 1 < imax) {

        imid = (int)((imin + imax) / 2);

        if (m_slow_solver[imid] <= axis_len) {

            imin = imid;

        } else {

            imax = imid;

        }

    }

    if (axis_len - m_slow_solver[imin] < m_slow_solver[imax] - axis_len) {

        return (double) imin / 10000;

    } else {

        return (double) imax / 10000;

    }

}

double Parameters::GetRadius(double axis_len) {

    if (axis_len <= 2 * pow(2,0.5)) {

        return axis_len / 2;

    } else {

       	double theta = m_fast_solver[HashAxisLength(axis_len)];
        return 4 * pow(M_PI, 0.5) / (sin(theta) - theta + 2 * M_PI);

    }

}

int Parameters::HashAxisLength(double axis_len) {

	return floor((axis_len - 2 * pow(2,0.5)) * 10000);

}

double Parameters::GetRandomGrowthRate() {

	return m_growth_dist[floor(R::runif(0, m_growth_dist.size()))];

}

double Parameters::GetMaxGrowth() {

	return *std::max_element(m_growth_dist.begin(), m_growth_dist.end());

}

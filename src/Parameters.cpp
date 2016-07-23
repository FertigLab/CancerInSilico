#include <cmath>
#include <Rcpp.h>
#include <algorithm>

#include "ExceptionHandling.h"
#include "Parameters.h"

void Parameters::InitializeRadiusSolver() {

	InitSlowSolver();
	InitFastSolver();

}

void Parameters::InitSlowSolver() {

    double theta;

    for (int i = 0; i <= 31400; ++i) {

        theta = (double) i / 10000;
        m_slow_solver.push_back(m_max_radius * pow(M_PI / (2 * M_PI - theta + sin(theta)), 0.5));

    }

}

void Parameters::InitFastSolver() {

	double theta;

	for (int i = 0; i <= 4142; ++i) {

		theta = GetThetaSlow(1 + (double) i / 10000);
		m_fast_solver.push_back(theta);

	}

}

double Parameters::GetThetaSlow(double rad) {

    int imin = 0, imax = 31400;
    int imid;

    while (imin + 1 < imax) {

        imid = (int)((imin + imax) / 2);

        if (m_slow_solver[imid] <= rad) {

            imin = imid;

        } else {

            imax = imid;

        }

    }

    if (rad - m_slow_solver[imin] < m_slow_solver[imax] - rad) {

        return (double) imin / 10000;

    } else {

        return (double) imax / 10000;

    }

}

double Parameters::GetTheta(double rad) {

	return m_fast_solver[HashRadius(rad)];

}

int Parameters::HashRadius(double rad) {

	return floor((rad - 1) * 10000);

}

double Parameters::GetRandomGrowthRate() {

	return m_growth_dist[floor(R::runif(0, m_growth_dist.size()))];

}

double Parameters::GetMaxGrowth() {

	return *std::max_element(m_growth_dist.begin(), m_growth_dist.end());

}

#include <cmath>
#include <Rcpp.h>
#include <algorithm>
#include <exception>

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
        num = pow(2 * M_PI, 0.5) * 2 * (1 + cos(theta / 2));
        denom = pow(sin(theta) - theta + 2 * M_PI, 0.5);
        m_slow_solver.push_back(num / denom);

    }

}

void Parameters::InitFastSolver() {

	double theta;

	for (int i = 0; i <= 11716; ++i) {

		theta = GetThetaSlow(2 * pow(2,0.5) + (double) i / 10000);
		m_fast_solver.push_back(theta);

	}

}

//m_slow_solver sorted from high to low
double Parameters::GetThetaSlow(double axis_len) {

    int imin = 0, imax = 31400;
    int imid;

    while (imin + 1 < imax) {

        imid = (int)((imin + imax) / 2);

        if (m_slow_solver[imid] >= axis_len) {

            imin = imid;

        } else {

            imax = imid;

        }

    }

    if (axis_len - m_slow_solver[imax] < m_slow_solver[imin] - axis_len) {

        return (double) imin / 10000;

    } else {

        return (double) imax / 10000;

    }

}

double Parameters::GetRadius(double axis_len) {

    if (axis_len < 2 * pow(2,0.5)) {

        throw std::invalid_argument("called deformation function with axis length less than max\n");

    } else {

       	double theta = m_fast_solver[HashAxisLength(axis_len)];
        return axis_len / (2 + 2 * cos(theta / 2));        

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

void Parameters::StoreGrowthDistribution(Rcpp::NumericVector growthRates) {

    for (unsigned int i = 0; i < growthRates.length(); ++i) {
    
        m_growth_dist.push_back(growthRates[i]);

    }

}

void Parameters::StoreDrugEffect(Rcpp::List drugEffect) {

    // Create DrugEffectMap
    DrugEffectMap drug_effect; 
    for (unsigned int i = 0; i < drugEffect.size(); ++i) {

        std::vector<double> dist = Rcpp::as< std::vector<double> >(drugEffect[i]);   

        double growthRate = dist[0];
        dist[0] = dist.back();
        dist.pop_back();     

        drug_effect.insert(std::pair<double, std::vector<double> >
                                (growthRate, dist));

    }

    m_drug_effect_map = drug_effect;

    // Insert comment here
    DrugEffectMap::iterator iter = m_drug_effect_map.begin();
    for (; iter != m_drug_effect_map.end(); ++iter) {

        m_drug_effect_indices.push_back(iter->first);

    }

    std::sort (m_drug_effect_indices.begin(), m_drug_effect_indices.end());

}


double Parameters::GetDrugEffect(double growthRate) {

    double index = *(std::lower_bound(m_drug_effect_indices.begin(),
                        m_drug_effect_indices.end(), growthRate));
    
    std::vector<double>& vec_ref = m_drug_effect_map.at(index);

    int size = vec_ref.size();
    return 1 - vec_ref[floor(R::runif(0,size))];
    
}
    

void Parameters::StoreCellTypes(Rcpp::List cellTypes) {

    for (unsigned int i = 0; i < cellTypes.size(); ++i) {
    
        m_cell_types.push_back(Rcpp::S4(cellTypes[i]));

    }

}

std::vector<Rcpp::S4> Parameters::GetCellTypes() {

    return m_cell_types;

}


void Parameters::StoreCellTypeDistribution(Rcpp::NumericVector cellTypeDist) {

    m_cell_type_dist = cellTypeDist;

}

int Parameters::GetRandomCellType() {

    Rcpp::IntegerVector cellTypeSample;

    R::rmultinom(1, m_cell_type_dist.begin(), 1, cellTypeSample.begin());

    for (unsigned int i = 0; i < cellTypeSample.length(); ++i) {
    
        if (cellTypeSample[i]) { // get the index
            return i;   
        }

    }

}
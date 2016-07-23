// [[Rcpp::depends(BH)]]

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <testthat.h>
#include <Rcpp.h>

class Parameters {

  private:

    int m_initial_num_cells;
    double m_max_translation, m_max_rotate, m_max_deform;
    double m_epsilon, m_delta;
	double m_max_radius;
	bool m_inherit_growth;
	double m_nG;

	std::vector<double> m_slow_solver;
	std::vector<double> m_fast_solver;
	std::vector<double> m_growth_dist;

    void InitializeRadiusSolver();
	void InitSlowSolver();
	void InitFastSolver();
	double GetThetaSlow(double);
	int HashRadius(double);

  public:

    Parameters(double max_rad) {
		m_max_radius = max_rad;
		InitializeRadiusSolver();
	}

    //Setters
    void SetMaxTranslation(double trans) { m_max_translation = trans;}
    void SetMaxRotate(double rot) { m_max_rotate = rot;}
    void SetMaxDeform(double def) { m_max_deform = def;}
    void SetResistanceEPSILON(double ep) { m_epsilon = ep;}
    void SetCompressionDELTA(double dt) { m_delta = dt;}
	void SetInheritGrowth(bool gr) { m_inherit_growth = gr;}
	void StoreGrowthDistribution(std::vector<double> gr) { m_growth_dist = gr;}
	void SetNG(double ng) { m_nG = ng;}

    //Getters
    double GetMaxTranslation() { return m_max_translation;}
    double GetMaxRotate() { return m_max_rotate;}
    double GetMaxDeform() { return m_max_deform;}
    double GetResistanceEPSILON() { return m_epsilon;}
    double GetCompressionDELTA() { return m_delta;}
	bool InheritGrowth() { return m_inherit_growth;}
	double GetNG() { return m_nG;}
	double GetMaxRadius() { return m_max_radius;}

	double GetRandomGrowthRate();
	double GetMaxGrowth();
    double GetTheta(double);

};

#endif

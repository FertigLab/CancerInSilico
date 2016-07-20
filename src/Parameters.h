// [[Rcpp::depends(BH)]]

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <testthat.h>
#include <Rcpp.h>

class Parameters {

  private:

    int m_initial_num_cells;
    double m_duration;
    double m_drug_prop;
    double m_max_migration, m_max_rotate, m_max_deform;
    double m_apoptosis_rate;
    double m_energy_constant;
    double m_epsilon, m_delta;
    double m_init_density;
	double m_max_radius;
	bool m_inherit_growth;

	std::vector<double> m_slow_solver;
	std::vector<double> m_fast_solver;
	std::vector<double> growth_rates;

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
    void SetInitialNumCells(int n) { m_initial_num_cells = n;}
    void SetDuration(double dur) { m_duration = dur;}
    void SetDrugProportion(double p) { m_drug_prop = p;}
    void SetMaxMigration(double mig) { m_max_migration = mig;}
    void SetMaxRotate(double rot) { m_max_rotate = rot;}
    void SetMaxDeform(double def) { m_max_deform = def;}
    void SetApoptosisRate(double apop) { m_apoptosis_rate = apop;}
    void SetEnergyConstant(double EC) { m_energy_constant = EC;}
    void SetResistanceEPSILON(double ep) { m_epsilon = ep;}
    void SetCompressionDELTA(double dt) { m_delta = dt;}
    void SetInitialDensity(double den) { m_init_density = den;}
	void SetInheritGrowth(bool gr) { m_inherit_growth = gr;}
	void StoreGrowthRates(std::vector<double> gr) { growth_rates = gr;}

    //Getters
    int GetInitialNumCells() { return m_initial_num_cells;}
    double GetDuration() { return m_duration;}
    double GetDrugProportion() { return m_drug_prop;}
    double GetMaxMigration() { return m_max_migration;}
    double GetMaxRotate() { return m_max_rotate;}
    double GetMaxDeform() { return m_max_deform;}
    double GetApoptosisRate() { return m_apoptosis_rate;}
    double GetEnergyConstant() { return m_energy_constant;}
    double GetResistanceEPSILON() { return m_epsilon;}
    double GetCompressionDELTA() { return m_delta;}
    double GetInitialDensity() { return m_init_density;}
	double GetMaxRadius() { return m_max_radius;}
	bool InheritGrowth() { return m_inherit_growth;}

	double GetRandomGrowthRate();
    double GetTheta(double);

};

#endif

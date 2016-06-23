#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <cmath>
#include <vector>

class Parameters {

private:

  int m_initial_num_cells;
  double m_duration;
  double m_min_radius, m_max_radius;
  double m_drug_prop;
  double m_mean_growth, m_var_growth;
  double m_max_migration, m_max_rotate, m_max_deform;
	double m_apoptosis_rate;
	double m_energy_constant;
	double m_epsilon, m_delta;
  double m_init_density;
  std::vector<double> m_radius_key;    

public:

  Parameters() {}

  //Setters
  void SetInitialNumCells(int n) {m_initial_num_cells = n;}
  void SetDuration(double dur) {m_duration = dur;}
  void SetMinRadius(double min_rad) {m_min_radius = min_rad;}
  void SetMaxRadius(double max_rad) {m_max_radius = max_rad;}
  void SetMeanGrowth(double mean) {m_mean_growth = mean;}
  void SetVarGrowth(double var) {m_var_growth = var;}
  void SetDrugProportion(double p) {m_drug_prop = p;}
  void SetMaxMigration(double mig) {m_max_migration = mig;}
  void SetMaxRotate(double rot) {m_max_rotate = rot;}
  void SetMaxDeform(double def) {m_max_deform = def;}
  void SetApoptosisRate(double apop) {m_apoptosis_rate = apop;}
  void SetEnergyConstant(double EC) {m_energy_constant = EC;}
  void SetResistanceEPSILON(double ep) {m_epsilon = ep;}
  void SetCompressionDELTA(double dt) {m_delta = dt;}
  void SetInitialDensity(double den) {m_init_density = den;}

  //Getters
  double GetMinRadius() {return m_min_radius;}
  double GetMaxRadius() {return m_max_radius;}
  double GetMeanGrowth() {return m_mean_growth;}
  double GetVarGrowth() {return m_var_growth;}
  int GetInitialNumCells() {return m_initial_num_cells;}
  double GetDuration() {return m_duration;}
  double GetDrugProportion() {return m_drug_prop;}
  double GetMaxMigration() {return m_max_migration;}
  double GetMaxRotate() {return m_max_rotate;}
  double GetMaxDeform() {return m_max_deform;}
  double GetApoptosisRate() {return m_apoptosis_rate;}
  double GetEnergyConstant() {return m_energy_constant;}
  double GetResistanceEPSILON() {return m_epsilon;}
  double GetCompressionDELTA() {return m_delta;}
  double GetNeighborRadius() {return m_max_radius + m_min_radius * 2;}
  double GetInitialDensity() {return m_init_density;}

  void InitializeRadiusSolver();
  double GetTheta(double);
  
};

#endif 

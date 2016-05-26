#ifndef CELL_POPULATION_HPP
#define CELL_POPULATION_HPP

#include <Rcpp.h>
#include <vector>
#include <cmath>

#include "Parameters.hpp"
#include "Cell.hpp"
#include "SpatialHash.hpp"

class CellPopulation {

private:

  Parameters* m_param;
  std::vector<Cell*> m_population;
  std::vector<std::vector<double> > m_population_record;

public:

  CellPopulation() {}
	~CellPopulation();
  CellPopulation(Parameters*, unsigned int, double);
  
  std::pair<double,double> GetRandomLocation(double);
  bool ValidCellPlacement(double,double);
  void OneTimeStep();
  void Update(SpatialHash&);
  void AttemptTrial(Cell*, SpatialHash&);
  bool AcceptTrial(double);
  double CalculateTotalInteraction(Cell*, SpatialHash&);
  double CalculateInteraction(Cell*,Cell*);
  void CheckMitosis();
	void RecordPopulation();
  void UpdateNeighbors(Cell*, SpatialHash&);
  Rcpp::NumericMatrix GetPopulationAsMatrix();
  void AddDrug();

};

#endif

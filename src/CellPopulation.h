#ifndef CELL_POPULATION_HPP
#define CELL_POPULATION_HPP

#include <Rcpp.h>
#include <vector>

#include "Point.h"
#include "Parameters.h"
#include "Cell.h"
#include "SpatialHash.h"

class CellPopulation {

friend class TestCellPopulation;

private:

	Parameters *m_param;

	SpatialHash<Cell> m_population;

	std::vector<std::vector<double> > m_population_record;

public:

	~CellPopulation();

	CellPopulation(Parameters *, unsigned int, double);

	Point GetRandomLocation(Cell*, double);
	bool ValidCellPlacement(Cell*);

	void OneTimeStep();
	void Update();

	void AttemptTrial(Cell *);
	bool AcceptTrial(double);
	double CalculateTotalInteraction(Cell *);
	double CalculateInteraction(Cell *, Cell *);

	void CheckMitosis(Cell*);
	void AddDrug();

	void RecordPopulation();
	Rcpp::NumericMatrix GetPopulationAsMatrix();

	int size();

};

#endif

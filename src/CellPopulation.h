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

    bool m_drug_added;

public:

	~CellPopulation();

	CellPopulation(Parameters*, unsigned int, double);

    void CreateBoundary();
    void CreateCells(int, double);
    void InitCellCycles();

	Point GetRandomLocation(Cell*, double);
	bool ValidCellPlacement(Cell*);

	void OneTimeStep();
	void Update();

	void AttemptTrial(Cell*);
    bool CheckForCellOverlap(Point, Cell*);
    bool CheckBoundary(Cell*);
    int CalculateNumberOfNeighbors(Cell*);
	double CalculateTotalInteraction(Cell*);
	double CalculateInteraction(Cell*, Cell*);
	bool AcceptTrial(double, double, Cell*);

	void CheckMitosis(Cell*);
	void SetGrowthRates();

    void AddDrug();

	void RecordPopulation();
	Rcpp::List GetPopulationAsList();

	int size();

};

#endif

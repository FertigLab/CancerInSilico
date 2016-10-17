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

    /* pointer to the parameters object */
	Parameters *m_param;

    /* data structure holding the population */
	SpatialHash<Cell> m_population;

    /* record of the cell information */
	std::vector<std::vector<double> > m_population_record;

    /* true once drug has been added to population */
    bool m_drug_added;

public:
    
    /* constructor; takes a Parameters object, initial size
       and density of  the population */
	CellPopulation(Parameters*, unsigned int, double);

    /* destructor */
	~CellPopulation();

    /****** Initialization Functions ******/

    /* create the cell boundary */
    void CreateBoundary(double);

    /* create the cells and seed them randomly inside the given radius */
    void CreateCells(int, double);

    /* find a random location inside the disk and move the cell there */
    void MoveToRandomLocation(Cell*, double);

    /* check if a cell is placed in a valid location */
    bool ValidCellPlacement(Cell*);

    /* set the initial growth rates of the cells */
    void SetGrowthRates();

    /* seed cells randomly throughout the cell cycle */
    void InitCellCycle();

    /* add the drug to the cell population */
    void AddDrug();

    /**************************************/

    /****** Functions for updating cells ******/

    /* take one model time step */
	void OneTimeStep();

    /* make one Monte Carlo step */
	void Update();

    /* attempt a growth/migration/rotation trial */
	void AttemptTrial(Cell*);

    /* check if cell overlaps any other cells near a point */
    bool CheckForCellOverlap(Point, Cell*);

    /* check if cell is within the boundary */
    bool CheckBoundary(Cell*);

    /* calculate how many neighbors the cell has */
    int CalculateNumberOfNeighbors(Cell*);

    /* calculate the interaction potential of this cell with all others */
	double CalculateTotalInteraction(Cell*);

    /* calculate the interaction potential between two cells */
	double CalculateInteraction(Cell*, Cell*);

    /* determine whether or not to accept the trail */
	bool AcceptTrial(double, double, Cell*);

    /* check if a cell is ready to divide */    
	void CheckMitosis(Cell*);

    /******************************************/

    /* record the current state of the population */    
	void RecordPopulation();

    /* return the cell population record as a R list */    
	Rcpp::List GetPopulationAsList();

    /* return the size of the population */
	int size();

};

#endif

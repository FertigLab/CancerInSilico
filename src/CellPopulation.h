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

    /* create a vector of 'default' state cells */
    std::vector<Cell*> CreateDefaultCells(unsigned int);

    /* seed cells randomly throughout the cell cycle */
    double InitCellCycle(std::vector<Cell*>);
    
    /* create the cell boundary */
    void CreateBoundary(double);

    /* create the cells and seed them randomly inside the given radius */
    void PlaceCells(std::vector<Cell*>, double);

    /* find a random location inside the disk and move the cell there */
    void MoveToRandomLocation(Cell*, double);

    /* check if a cell is placed in a valid location */
    bool ValidCellPlacement(Cell*, double);

    /* set the initial growth rates of the cells */
    void SetGrowthRates();

    /* add the drug to the cell population */
    void AddDrug();

    /**************************************/

    /****** Functions for updating cells ******/

    /* take one model time step */
	void OneTimeStep();

    /* take one Monte Carlo step */
	void Update();

    /* attempt a growth/migration/rotation trial */
	void AttemptTrial(Cell*);

    /* calculate how many neighbors the cell has */
    int CalculateNumberOfNeighbors(Cell*);

    /* check if cell overlaps any other cells near a point */
    bool CheckForCellOverlap(Point, Cell*);

    /* check if cell is within the boundary */
    bool CheckBoundary(Cell*, double);

    /* determine whether or not to accept the trial */
	bool AcceptTrial(double, double, Cell*);

    /* calculate the interaction potential of this cell with all others */
	double CalculateTotalInteraction(Cell*);

    /* calculate the interaction potential between two cells */
	double CalculateInteraction(Cell*, Cell*);

    /* check if a cell is ready to divide */    
	void CheckMitosis(Cell*);

    /******************************************/

    /* record the current state of the cell population */    
	void RecordPopulation();

    /* return the cell population record as an R list */    
	Rcpp::List GetPopulationAsList();

    /* return the size of the cell population */
    int size();

};

#endif

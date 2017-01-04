#ifndef CELL_POPULATION_HPP
#define CELL_POPULATION_HPP

#include <Rcpp.h>
#include <vector>

#include "Point.h"
#include "Parameters.h"
#include "Cell.h"
#include "SpatialHash.h"

class CellPopulation
{

friend class TestCellPopulation;

private:

    /* pointer to the parameters object */
	Parameters* mParameters;

    /* data structure holding the cell population */
	SpatialHash<Cell> mPopulation;

    /* record of the cell information */
	std::vector< std::vector<double> > mPopulationRecord;

    /* create a vector of 'default' state cells */
    std::vector<Cell*> createDefaultCells(unsigned int);

    /* calc boundary, return seeding radius - different if no boundary */
    double calculateBoundary(std::vector<Cell*>&);

    /* create the cells and seed them randomly inside the given radius */
    void placeCells(std::vector<Cell*>&, double);

public:
    
    /* constructor; takes a Parameters object, initial size
       and density of  the population */
	CellPopulation(Parameters*);

    /* destructor */
	~CellPopulation();

    /* find a random location inside the disk and move the cell there */
    void placeRandomly(Cell*);

    /* check if cell is in valid location */
    bool validCellLocation(Cell*);

    /* calculate how many neighbors the cell has */
    int numberOfNeighbors(Cell*);

    /* record the current state of the cell population */    
	void recordPopulation();

    /* return the cell population record as an R list */    
	Rcpp::List populationAsList();

    /* return the size of the cell population */
    int size();
};

#endif

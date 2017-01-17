#include <cmath>
#include <Rcpp.h>

#include "CellPopulation.h"

/****** CONSTRUCTOR/DESTRUCTOR AND INITIALIZATION FUNCTIONS ******/

/* constructor; takes a Parameters object, initial size and
   density of the  population */
CellPopulation::CellPopulation(Parameters *par, unsigned int size, double density) {

    /* store the parameters */
    m_param = par;

    /* create a temporary vector of 'default' state cells */
    std::vector<Cell*> cells = CreateDefaultCells(size);

    /* seed cells throughout the cell cycle */
    double cell_area = InitCellCycle(cells);

    /* calculate the initial seeding radius */
    double disk_radius = pow(cell_area / (M_PI * density), 0.5);
   
    /* create the cell boundary: cells cannot go outside the boundary */
    CreateBoundary(disk_radius);

    /* initalize the data structure to hold the cells */
    m_population = SpatialHash<Cell>(1.0);
  
    /* place the cells and seed them randomly throughout the disk;
       cells are inserted in m_population during this step */
    PlaceCells(cells, disk_radius);

    /* intialize to false: haven't added the drug yet */
    m_drug_added = false;

}

/* deconstructor */
CellPopulation::~CellPopulation() {
    
    /* iterate through entire cell population */
    SpatialHash<Cell>::full_iterator iter = m_population.begin();
    for (; iter != m_population.end(); ++iter) {

        /* free memory of cell */
        delete &iter;

    }

}

/* create a temporary vector of 'default' state cells */
std::vector<Cell*> CellPopulation::CreateDefaultCells(unsigned int num) {

    /* create the return vector */
    std::vector<Cell*> ret_vec;

    /* create 'num' cells */
    for (unsigned int i = 0; i < num; ++i) {

        /* use default constructor at point (0,0) */
        ret_vec.push_back(new Cell(Point(0,0), m_param));

    }

    /* return cells */
    return ret_vec;

}

/* seed cells randomly throughout the cell cycle, return total area
   occupied by the cells */
double CellPopulation::InitCellCycle(std::vector<Cell*> cells) {

    /* declare variable for later */
    double unif, total_area = 0.0;

    /* get iterator */
    std::vector<Cell*>::iterator it = cells.begin();

    /* iterate through each cell */
    for (; it != cells.end(); ++it) {

        /* random number in (0,1) for probability calculations */        
        unif = R::runif(0,1);

        /* probability of being seeded in interphase */
        if (m_param->syncCycles() || unif < 0.75) { //interphase

            /* set random radius and resulting axis length */
            (*it)->SetRadius(R::runif(1,m_param->maxRadius()));
            (*it)->SetAxisLength((*it)->GetRadius() * 2);
    
        /* otherwise seed in mitosis */
        } else { //mitosis

            /* put the cell in a random point in mitosis */
            (*it)->EnterRandomPointOfMitosis();

        }

        /* add area of cell to total */
        total_area += (*it)->GetArea();

    }

    /* return the total area the cells occupy */
    return total_area;

}


/* create the cell boundary */
void CellPopulation::CreateBoundary(double radius) {

    /* if boundary is zero (no boundary) */
    if (m_param->boundary() == 0.0) {

        /* set boundary value to maximum allowed */
        m_param->setBoundary(std::numeric_limits<double>::max());

    /* otherwise check if boundary is too small
       (must be bigger than initial seeding radius) */
    } else if (m_param->boundary() < radius) {

        /* set boundary to the minimum value */
        m_param->setBoundary(radius);

    }

}

/* create all the cells */
void CellPopulation::PlaceCells(std::vector<Cell*> cells, double radius) {

    /* get cell iterator */
    std::vector<Cell*>::iterator it = cells.begin();

    /* create 'num' amount of cells */
    for (; it != cells.end(); ++it) {

        /* get random location inside the radius and move the cell there */
        MoveToRandomLocation(*it, radius);

        /* add cell to the cell population */
        m_population.Insert((*it)->GetCoord(), *it);

        /* check if the user wants to cancel the simulation
           (done from R console) */
        Rcpp::checkUserInterrupt();

    }

}

/* find a random location inside the disk and move the cell there */
void CellPopulation::MoveToRandomLocation(Cell* cell, double rad) {

    /* delcare variables used later */
    double dist, ang, x, y;

    do {

        /* calculate a random distance and angle */    
        dist = R::runif(0, 1);
        ang = R::runif(0, 2 * M_PI);

        /* calculate point coordinates */
        x = rad * pow(dist, 0.5) * cos(ang);
        y = rad * pow(dist, 0.5) * sin(ang);

        /* move the cell to the new coordinates */
        cell->SetCoord(Point(x,y));

        /* check if the user wants to cancel the simulation
           (done from R console) */
        Rcpp::checkUserInterrupt();

    /* keep trying points until the cell finds a valid location */
    } while (!ValidCellPlacement(cell, rad));

}

/* check if a cell is placed in a valid location */
bool CellPopulation::ValidCellPlacement(Cell* cell, double rad) {

    /* check if cell exceeds seeding radius */
    if (CheckBoundary(cell, rad)) {

        return false;

    }

    /* get iterator for entire population */
    SpatialHash<Cell>::full_iterator iter = m_population.begin();

    /* iterate through whole population */
    for (; iter != m_population.end(); ++iter) {

        /* check if this cell overlaps */
        if ((*iter).CellDistance(*cell) < 0) {

            /* return false if overlap is found */
            return false;

        }

    }

    /* if no overlap found, return true */
    return true;

}

/* add the drug to the cell population */
void CellPopulation::AddDrug() {

    /* get a iterator to the entire population */    
    SpatialHash<Cell>::full_iterator iter = m_population.begin();
    
    /* iterate through every cell */
    for (; iter != m_population.end(); ++iter) {

        /* get the current growth rate */
        double gr = (*iter).GetGrowth();
    
        /* adjust growth rate based on drug effect */
        (*iter).SetGrowth(gr * m_param->GetDrugEffect(gr));

    }

    /* record that the drug was added */
    m_drug_added = true;

}

/*****************************************************************/

/* take one model time step */
void CellPopulation::OneTimeStep() {

    /* store size of cell population */
    int sz = size();

    /* take 'sz' number of Monte Carlo steps */
    for (int i = 0; i < sz; ++i) {

        /* take a single Monte Carlo step */
        Update();

    }

}

/* take one Monte Carlo step */
void CellPopulation::Update() {

    /* get a random cell */    
    Cell* rand_cell = m_population.GetRandomValue();

    /* attempt a trial */
    AttemptTrial(rand_cell);

    /* check if the cell needs to divide */
    CheckMitosis(rand_cell);

}

/* attempt a growth/migration/rotation trial */
void CellPopulation::AttemptTrial(Cell *cell) {

    /* get interaction potential before the trial */
    double interaction = CalculateTotalInteraction(cell);

    /* get the number of neighbors around the cell */
    int num_neighbors = CalculateNumberOfNeighbors(cell);

    /* store the original cell */    
    Cell orig = *cell;

    /* attempt a trial */
    bool growth = cell->DoTrial();

    /* check if cell now overlaps any other */
    bool overlap = CheckForCellOverlap(orig.GetCoord(), cell);

    /* if the cell overlaps or crossed the boundary */
    if (overlap || CheckBoundary(cell, m_param->boundary())) {

        /* reject the trial, revert cell to original state */
        *cell = orig;

    /* if the trial wasn't a growth trial (automatic accept) */
    } else if (!growth) {

        /* update data structure with cells new position */
        m_population.Update(orig.GetCoord(), cell->GetCoord());    
        
        /* determine if trial is accepted */
        if (!AcceptTrial(interaction, num_neighbors, cell)) {

            /* if rejected, revert data struture and cell */
            m_population.Update(cell->GetCoord(), orig.GetCoord());    
            *cell = orig;

        }

    }

}

/* calculate how many neighbors the cell has */
int CellPopulation::CalculateNumberOfNeighbors(Cell *cell) {

    /* initalize count */
    int neighbors = 0;

    /* serach radius, all cells in this radius could have a
       non-zero interaction with the current cell (defn of neighbor) */
    double radius = m_param->delta() + 1;
    
    /* iterator in search radius */
    SpatialHash<Cell>::circular_iterator iter
        = m_population.begin(cell->GetCoord(), radius);

    /* step through iterator */    
    for (; iter != m_population.end(cell->GetCoord(), radius); ++iter) {

        /* check if cell is indeed inside this radius */
        if (cell->CellDistance(*iter) <= m_param->delta()) {

            /* increment the neighbor count */
            neighbors++;

        }

    }

    /* return the number of neighbors */
    return neighbors;

}

/* check if cell overlaps any other cells near a point */
bool CellPopulation::CheckForCellOverlap(Point center, Cell* cell) {

    /* calculate maximum search radius for overlapping cells */
    double max_search = std::max(m_param->maxTranslation(),
        m_param->maxRadius());

    /* get iterator for cells in search radius */
    SpatialHash<Cell>::circular_iterator iter
        = m_population.begin(center, max_search);

    /* iterate through cells */
    for (; iter != m_population.end(center, max_search); ++iter) {

        /* return true if cells overlap */
        if ((*iter).CellDistance(*cell) < 0) {

            return true;

        }

    }

    /* if no overlap found, return false */
    return false;

}

/* check if cell is within the boundary */
bool CellPopulation::CheckBoundary(Cell* cell, double bound) {

    /* get x,y offset of each part of the cell */
    double x_dist = (0.5 * cell->GetAxisLength() - cell->GetRadius())
        * cos(cell->GetAxisAngle());
    double y_dist = (0.5 * cell->GetAxisLength() - cell->GetRadius())
        * sin(cell->GetAxisAngle());

    /* calculate the center of each part of the cell */
    Point center_1 = Point(cell->GetCoord().x + x_dist,
        cell->GetCoord().y + y_dist);
    Point center_2 = Point(cell->GetCoord().x - x_dist,
        cell->GetCoord().y - y_dist);

    /* return true if cell is farther from center than the boundary line */
    return (center_1.dist(Point(0,0)) + cell->GetRadius() > bound
        || center_2.dist(Point(0,0)) + cell->GetRadius() > bound);

}

/* determine whether or not to accept the trial */
bool CellPopulation::AcceptTrial(double prev_interaction,
    double prev_num_neighbors, Cell* cell) {

    /* calculate the difference in interaction potentials */
    double delta_interaction = CalculateTotalInteraction(cell)
        - prev_interaction;

    /* see if cell lost a neighbor */
    if (CalculateNumberOfNeighbors(cell) < prev_num_neighbors) {

        /* automatic rejection */
        return false;

    /* if interaction potential decreased */
    } else if (delta_interaction <= 0.0) {

        /* automatic acceptance */
        return true;

    /* otherwise accept with certain probability */
    } else {

        /* get random number in (0,1) for probability calculation */
        double unif = R::runif(0, 1);

        /* calculate probability of acceptance */
        double prob = exp(-1 * delta_interaction);

        /* return true with probability 'prob' */
        return unif < prob;

    }

}

/* calculate the interaction potential of this cell with all others */
double CellPopulation::CalculateTotalInteraction(Cell *cell) {

    /* initalize total interaction to zero */
    double sum = 0.0;

    /* get search radius for cells that could have positive interactions */
    double rad = m_param->delta() + 1;

    /* get iterator for cells in this radius */    
    SpatialHash<Cell>::circular_iterator iter
        = m_population.begin(cell->GetCoord(), rad);

    /* iterate through cells */
    for (; iter != m_population.end(cell->GetCoord(), rad); ++iter) {

        /* add the individual interactions together */
        sum += CalculateInteraction(&iter, cell);

    }

    /* return total interaction */
    return sum;

}

/* calculate the interaction potential between two cells */
double CellPopulation::CalculateInteraction(Cell* a, Cell* b) {

    /* get distance between the two cells */
    double dist = a->CellDistance(*b);

    /* if distance is too far, no interaction happens */
    if (dist > m_param->delta()) {

        return 0.0;
    
    /* otherwise calculate the interaction */
    } else {

        /* get one component of the interaction */
        double part = 2 * dist / m_param->delta() - 1;

        /* calculate interaction */
        return m_param->epsilon() * (pow(part,2) - 1);

    }

}

/* check if a cell is ready to divide */
void CellPopulation::CheckMitosis(Cell* cell) {

    /* if cell is ready to divide */
    if (cell->ReadyToDivide()) {

        /* get coordinates of parent cell */
        Point old_key = cell->GetCoord();

        /* create daughter cell */
        Cell* daughter_cell = new Cell(cell->Divide(m_drug_added));

        /* add daughter cell to population */
        m_population.Insert(daughter_cell->GetCoord(), daughter_cell);

        /* update parents coordinates to new position */
        m_population.Update(old_key, cell->GetCoord());

    }

}

/* record the current state of the cell population */
void CellPopulation::RecordPopulation() {

    /* the vector to hold the current population */
    std::vector<double> current_pop;
    
    /* cell population iterator */
    SpatialHash<Cell>::full_iterator iter = m_population.begin();
  
    /* loop through each cell */  
    for (; iter != m_population.end(); ++iter) {

        /* store cell information */
        current_pop.push_back((*iter).GetCoord().x);
        current_pop.push_back((*iter).GetCoord().y);
        current_pop.push_back((*iter).GetRadius());
        current_pop.push_back((*iter).GetAxisLength());
        current_pop.push_back((*iter).GetAxisAngle());
        current_pop.push_back((*iter).GetGrowth());
        current_pop.push_back((*iter).GetType());

    }

    /* add current population to record */    
    m_population_record.push_back(current_pop);

}

/* return the cell population record as an R list */
Rcpp::List CellPopulation::GetPopulationAsList() {

    /* convert population record to R list */
    return Rcpp::wrap(m_population_record);

}

/* return the size of the cell population */
int CellPopulation::size() {

    /* return size */
    return m_population.size();

}

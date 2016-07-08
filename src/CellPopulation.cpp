#include <cmath>
#include <Rcpp.h>

#include "CellPopulation.h"

CellPopulation::CellPopulation(Parameters *par, unsigned int size, double density) {

    m_param = par;
    m_population = SpatialHash<Cell>(1.0);
    double disk_radius = pow(size / density, 0.5);
    Point new_loc;
	Cell* temp;

    for (unsigned int i = 0; i < size; i++) {

        new_loc = GetRandomLocation(disk_radius);
		temp = new Cell(new_loc, m_param);
        m_population.Insert(temp->GetCoord(), temp);

    }

}

CellPopulation::~CellPopulation() {
	
	SpatialHash<Cell>::full_iterator iter = m_population.begin();
    for (; iter != m_population.end(); ++iter) {

        delete &iter;

    }

}

Point CellPopulation::GetRandomLocation(double rad) {

    double dist, ang, x, y;

    do {

        dist = R::runif(0, 1);
        ang = R::runif(0, 2 * M_PI);
        x = rad * pow(dist, 0.5) * cos(ang);
        y = rad * pow(dist, 0.5) * sin(ang);

    } while (!ValidCellPlacement(Point(x,y)));

    return Point(x, y);

}

bool CellPopulation::ValidCellPlacement(Point pt) {

    SpatialHash<Cell>::full_iterator iter = m_population.begin();

	for (; iter != m_population.end(); ++iter) {

        if ((*iter).CellDistance(Cell(pt, m_param)) < 0) {

            return false;

        }

    }

    return true;

}

void CellPopulation::OneTimeStep() {

    int sz = m_population.size();

    for (int i = 0; i < sz; ++i) {

        Update();

    }

    RecordPopulation();

}

void CellPopulation::Update() {

    Cell* rand_cell = m_population.GetRandomValue();
    AttemptTrial(rand_cell);
	CheckMitosis(rand_cell);

}

void CellPopulation::CheckMitosis(Cell* cell) {

	if (cell->ReadyToDivide()) {
	
		Point old_key = cell->GetCoord();
		Cell* daughter_cell = new Cell(cell->Divide());
		m_population.Insert(daughter_cell->GetCoord(), daughter_cell);
		m_population.Update(old_key, cell->GetCoord());

	}

}

void CellPopulation::AttemptTrial(Cell *cell) {

    double pre_interaction = CalculateTotalInteraction(cell);
    Cell orig = *cell;
    cell->DoTrial();

	double max_search = std::max(m_param->GetMaxMigration(), m_param->GetMaxRadius());

	SpatialHash<Cell>::circular_iterator iter
		= m_population.begin(orig.GetCoord(), max_search);


	for (; iter != m_population.end(orig.GetCoord(), max_search); ++iter) {

		if ((*iter).CellDistance(*cell) < 0) {

			*cell = orig;

		}

	}

    m_population.Update(orig.GetCoord(), cell->GetCoord());

    double post_interaction = CalculateTotalInteraction(cell);

    if (post_interaction == std::numeric_limits<double>::max()
            || !AcceptTrial(post_interaction - pre_interaction)) {

        m_population.Update(cell->GetCoord(), orig.GetCoord());
        *cell = orig;

    }

}

bool CellPopulation::AcceptTrial(double delta_interaction) {

    if (delta_interaction <= 0.0) {

        return true;

    } else {

        double unif = R::runif(0, 1);
        double prob = exp(-1 * delta_interaction / m_param->GetEnergyConstant());
        return unif < prob;

    }

}

double CellPopulation::CalculateTotalInteraction(Cell *cell) {

    double inter, sum = 0.0;

	SpatialHash<Cell>::circular_iterator iter
		= m_population.begin(cell->GetCoord(), m_param->GetCompressionDELTA() + 1);

	for (; iter != m_population.end(cell->GetCoord(), m_param->GetCompressionDELTA() + 1); ++iter) {

        if (&iter != cell) {

            inter = CalculateInteraction(&iter, cell);

            if (inter == std::numeric_limits<double>::max()) {

                return inter;

            }

            sum += inter;

        }

    }

    return sum;

}

double CellPopulation::CalculateInteraction(Cell* a, Cell* b) {

    double dist = a->CellDistance(*b);

    if (dist > m_param->GetCompressionDELTA()) {

        return 0.0;

    } else if (dist < 0) { //should never be called

        return std::numeric_limits<double>::max();

    } else {

        double part = pow(2 * dist / m_param->GetCompressionDELTA(), 2);
        return m_param->GetResistanceEPSILON() * (part - 1);

    }

}


void CellPopulation::AddDrug() {

    double rand;
	SpatialHash<Cell>::full_iterator iter = m_population.begin();

	for (; iter != m_population.end(); ++iter) {

        rand = R::rnorm(m_param->GetMeanGrowth(), m_param->GetVarGrowth());
        (*iter).SetGrowth(std::max(rand, 0.02));

    }

}

int CellPopulation::size() {

    return m_population.size();

}

void CellPopulation::RecordPopulation() {

    std::vector<double> current_pop;
	SpatialHash<Cell>::full_iterator iter = m_population.begin();
    Cell *temp;

	for (; iter != m_population.end(); ++iter) {

        current_pop.push_back((*iter).GetCoord().x);
        current_pop.push_back((*iter).GetCoord().y);
        current_pop.push_back((*iter).GetRadius());
        current_pop.push_back((*iter).GetAxisLength());
        current_pop.push_back((*iter).GetAxisAngle());
        current_pop.push_back((*iter).GetGrowth());

    }

    m_population_record.push_back(current_pop);

}

Rcpp::NumericMatrix CellPopulation::GetPopulationAsMatrix() {

    int num_rows = m_population_record.size();
    unsigned int num_cols = 0;

    for (int i = 0; i < num_rows; i++) {
 
       if (m_population_record[i].size() > num_cols) {

            num_cols = m_population_record[i].size();

        }

    }

    Rcpp::NumericMatrix ret_val(num_rows, num_cols);

    for (int i = 0; i < num_rows; i++) {

        for (unsigned int j = 0; j < m_population_record[i].size(); j++) {

            ret_val(i, j) = m_population_record[i][j];

        }

    }

    return ret_val;

}


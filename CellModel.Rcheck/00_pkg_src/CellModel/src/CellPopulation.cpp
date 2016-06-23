#include <cmath>
#include <cstdlib>
#include <Rcpp.h>

#include "CellPopulation.hpp"
#include "SpatialHash.hpp"
#include "Cell.hpp"

CellPopulation::CellPopulation(Parameters* par, unsigned int size, double density) {

  m_param = par;
  m_param->InitializeRadiusSolver();
  m_population = SpatialHash(m_param->GetMinRadius());
  
  double disk_radius = pow(size * pow(m_param->GetMinRadius(),2) / density, 0.5);
  std::pair<double,double> new_loc;

  for (unsigned int i = 0; i < size; i++) {

    new_loc = GetRandomLocation(disk_radius);
    m_population.Insert(new Cell(new_loc, m_param));   
    
  }

}

CellPopulation::~CellPopulation() {

/*  for (unsigned int i = 0; i < m_population.size(); i++) {
    
    delete m_population[i];
  
  }*/

}

std::pair<double,double> CellPopulation::GetRandomLocation(double rad) {

  double dist,ang,x,y;

  do {
    
    dist = R::runif(0,1);  
    ang = R::runif(0,2 * M_PI);  
    x = rad * pow(dist,0.5) * cos(ang);
    y = rad * pow(dist,0.5) * sin(ang);

  } while (!ValidCellPlacement(x,y));

  return std::make_pair(x,y);

}

bool CellPopulation::ValidCellPlacement(double x, double y) {

  SpatialIterator iter = m_population.getFullIterator();

  while (iter.Next()) {

    Cell temp (std::make_pair(x,y), m_param);
    if (iter.getCell()->CellDistance(temp) < 0) {

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

  Cell* rand_cell = m_population.GetRandomCell();

  AttemptTrial(rand_cell);
  CheckMitosis();

}


void CellPopulation::AttemptTrial(Cell* cell) {
  
  double pre_interaction = CalculateTotalInteraction(cell);  
  
  Cell orig = *cell;

  try {
  
    m_population.Update(orig, cell->DoTrial());

  } catch (std::exception& e) {

    *cell = orig;
    m_population.AddKey(cell);

  }

  double post_interaction = CalculateTotalInteraction(cell);
  
  if (post_interaction == std::numeric_limits<double>::max()
        || !AcceptTrial(post_interaction - pre_interaction)) {
    
    m_population.Update(*cell, orig);
    *cell = orig;
    
  }

}

bool CellPopulation::AcceptTrial(double delta_interaction) {
  
  if (delta_interaction <= 0.0) {

    return true;

  } else {

    double unif = R::runif(0,1);
    double prob = exp(-1 * delta_interaction / m_param->GetEnergyConstant());
    return unif < prob;

	}

}

double CellPopulation::CalculateTotalInteraction(Cell* cell) {
  
  double inter, sum = 0.0;  
  SpatialIterator iter = m_population.getCircularIterator(cell, m_param->GetCompressionDELTA() + 1);
  
  while (iter.Next()) {

    if (*(iter.getCell()) != *cell) {

      inter = CalculateInteraction(iter.getCell(), cell);
			if (inter == std::numeric_limits<double>::max()) {return inter;}
      sum += inter;

    }

  }
  
  return sum;

}

double CellPopulation::CalculateInteraction(Cell* a, Cell* b) {

  double dist = a->CellDistance(*b);

  if (dist > m_param->GetCompressionDELTA()) {  

    return 0.0;

  } else if (dist < 0) {

		return std::numeric_limits<double>::max();

  } else {

    double part = pow(2 * dist / m_param->GetCompressionDELTA(),2);
    return m_param->GetResistanceEPSILON() * (part - 1);

  }

}


void CellPopulation::CheckMitosis() {

  std::vector<Cell*> dividing_cells;
  SpatialIterator iter = m_population.getFullIterator();  

  while (iter.Next()) {

    if (iter.getCell()->ReadyToDivide()) {

			dividing_cells.push_back(iter.getCell());

		}

	}

  for (unsigned int i = 0; i < dividing_cells.size(); ++i) {

    Cell temp = *dividing_cells[i];
    m_population.Insert(dividing_cells[i]->Divide());
    m_population.Update(temp, *dividing_cells[i]);

  }

}

void CellPopulation::RecordPopulation() {

  std::vector<double> current_pop;
  SpatialIterator iter = m_population.getFullIterator();

  Cell* temp;

  while (iter.Next()) {

    temp = iter.getCell();
    current_pop.push_back(temp->GetCoord().first);		
    current_pop.push_back(temp->GetCoord().second);
    current_pop.push_back(temp->GetRadius());
    current_pop.push_back(temp->GetAxisLength());
    current_pop.push_back(temp->GetAxisAngle());
    current_pop.push_back(temp->GetGrowth());
    
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

      ret_val(i,j) = m_population_record[i][j];

    }

  }

  return ret_val;

}

void CellPopulation::AddDrug() {
  
  double rand;
  SpatialIterator iter = m_population.getFullIterator();

  while (iter.Next()) {

    rand = R::rnorm(m_param->GetMeanGrowth(), m_param->GetVarGrowth());
    iter.getCell()->SetGrowth(std::max(rand,0.02));

  }

}

int CellPopulation::size() {

  return m_population.size();

}






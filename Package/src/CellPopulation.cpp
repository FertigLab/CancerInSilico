#include <cmath>
#include <cstdlib>

#include <Rcpp.h>
#include "CellPopulation.hpp"
#include "Cell.hpp"

#include <iostream>

CellPopulation::CellPopulation(Parameters* par, unsigned int size, double density) {

  m_param = par;
  m_param->InitializeRadiusSolver();
  
  double disk_radius = pow(size * pow(m_param->GetMinRadius(),2) / density, 0.5);
  Cell* temp;

  for (int i = 0; i < size; i++) {

    std::pair<double,double> new_loc = GetRandomLocation(disk_radius);
    m_population.push_back(new Cell(new_loc, m_param));   
    
  }
 
}

CellPopulation::~CellPopulation() {

  std::vector<Cell*>::iterator iter = m_population.begin();

  for (; iter != m_population.end(); ++iter) {
    
    delete *iter;

  }


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

  std::vector<Cell*>::iterator iter = m_population.begin();

  for (; iter != m_population.end(); ++iter) {

    Cell temp (std::make_pair(x,y), m_param);
    if ((*iter)->CellDistance(temp) < 0) {

      return false;

    }

  }

  return true;

}

void CellPopulation::OneTimeStep() {

  int sz = m_population.size();
  SpatialHash hash_map = SpatialHash(m_population);  

  for (int i = 0; i < sz; i++) {

    Update(hash_map);

  }
  
	RecordPopulation();

}

void CellPopulation::Update(SpatialHash& hash){

  Cell* rand_cell = m_population[rand() % m_population.size()];
  UpdateNeighbors(rand_cell, hash); 
  AttemptTrial(rand_cell, hash);
  CheckMitosis();

}

void CellPopulation::UpdateNeighbors(Cell* cell, SpatialHash& hash) {

  double dist, min_dist = std::numeric_limits<double>::max();
  BucketIterator iter = hash.getCircularIterator(cell, 4);
  
  while (iter.Next()) {
    
    if (*(iter.getCell()) != *cell) {    
      
      dist = cell->CellDistance(*(iter.getCell()));

      if (dist < min_dist) {

        min_dist = dist;

      }

    }

  }

  cell->SetMinCellDist(min_dist);

}

void CellPopulation::AttemptTrial(Cell* cell, SpatialHash& hash) {
  
  Cell temp = *cell;
 
  double pre_interaction = CalculateTotalInteraction(cell, hash);  
  cell->DoTrial();
	double post_interaction = CalculateTotalInteraction(cell, hash);
  
  if (post_interaction == std::numeric_limits<double>::max()
        || !AcceptTrial(post_interaction - pre_interaction)) {

    *cell = temp;

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

double CellPopulation::CalculateTotalInteraction(Cell* cell, SpatialHash& hash) {
  
  double inter, sum = 0.0;  
  BucketIterator iter = hash.getCircularIterator(cell, m_param->GetCompressionDELTA() + 1);
  
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

  int sz = m_population.size();

  for (unsigned int i = 0; i < sz; ++i) {

    if (m_population[i]->ReadyToDivide()) {

			m_population.push_back(m_population[i]->Divide());

		}

	}

}

void CellPopulation::RecordPopulation() {

  std::vector<double> current_pop;
  std::vector<Cell*>::iterator iter = m_population.begin();

  for (; iter != m_population.end(); ++iter) {

		std::pair<double, double> crd = (*iter)->GetCoord();
    current_pop.push_back(crd.first);		
    current_pop.push_back(crd.second);
    current_pop.push_back((*iter)->GetRadius());
    current_pop.push_back((*iter)->GetAxisLength());
    current_pop.push_back((*iter)->GetAxisAngle());
    current_pop.push_back((*iter)->GetGrowth());
    
  }
  
  m_population_record.push_back(current_pop);

}

Rcpp::NumericMatrix CellPopulation::GetPopulationAsMatrix() {

  int num_rows = m_population_record.size();
  int num_cols = 0;

  for (int i = 0; i < num_rows; i++) {
    
    if (m_population_record[i].size() > num_cols) {

      num_cols = m_population_record[i].size();

    }

  }
  
  Rcpp::NumericMatrix ret_val(num_rows, num_cols);

  for (int i = 0; i < num_rows; i++) {

    for (int j = 0; j < m_population_record[i].size(); j++) {

      ret_val(i,j) = m_population_record[i][j];

    }

  }

  return ret_val;

}

void CellPopulation::AddDrug() {
  
  double rand;
  std::vector<Cell*>::iterator iter = m_population.begin();

  for (; iter != m_population.end(); ++iter) {

    rand = R::rnorm(m_param->GetMeanGrowth(), m_param->GetVarGrowth());
    (*iter)->SetGrowth(std::max(rand,0.02));

  }

}







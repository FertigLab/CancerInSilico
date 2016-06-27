#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Rcpp.h>

#include "Cell.h"

//used only for the initial population of cells
Cell::Cell(std::pair<double, double> coor, Parameters* par) {

  m_param = par;
  m_coordinates = coor;
  m_in_mitosis = false;
  m_ready_to_divide = false;  
  m_axis = std::make_pair(0.0, 0.0);
  m_radius = R::runif(m_param->GetMinRadius(), m_param->GetMaxRadius());
  m_growth_rate = m_param->GetMeanGrowth();

}

//used only for daughter cells
Cell::Cell(std::pair<double, double> coor, Parameters* par, double gr_rate) {

  m_param = par;
  m_coordinates = coor;
  m_in_mitosis = false;
  m_ready_to_divide = false;  
  m_axis = std::make_pair(0.0, 0.0);
  m_radius = m_param->GetMinRadius();
  m_growth_rate = gr_rate;

}

Cell& Cell::DoTrial() {

  double unif = R::runif(0,1);

  if (!m_in_mitosis) { //Interphase

    if (unif <= 0.5) {
      Migration();
    } else {
      Growth();
    }

  } else { //Mitosis

    if (unif <= 0.333) {
      Migration();
    } else if (unif <= 0.666) {
      Rotation();
    } else if (!m_ready_to_divide) {
      Deformation();
    }

  }

  return *this;

}
 
void Cell::Migration() {
  
  double length = m_param->GetMaxMigration() * pow(R::runif(0,1),0.5);
  double direction = R::runif(0,2 * M_PI);

  m_coordinates.first += length * cos(direction);
  m_coordinates.second += length * sin(direction);

}

void Cell::Growth() {
  
  double max_growth = 0.01 + m_param->GetMaxRadius() - m_radius;
  double growth = R::runif(0,std::min(m_growth_rate,max_growth));
  m_radius += growth;
  
  if (m_radius >= m_param->GetMaxRadius()) {

    m_radius = m_param->GetMaxRadius();  
    m_axis.first = 2 * m_radius;  
    m_in_mitosis = true;

  }

}

void Cell::Rotation() {

  double rotate = R::runif(-m_param->GetMaxRotate(), m_param->GetMaxRotate());
  m_axis.second += rotate;

}


void Cell::Deformation() {

  double max_deform = std::min(m_radius - m_param->GetMinRadius() + 0.01, m_param->GetMaxDeform());
  double deform = R::runif(0,max_deform);

  m_radius = std::max(m_radius - deform, m_param->GetMinRadius());
  m_axis.first = 2 * m_radius * (1 + cos(m_param->GetTheta(m_radius) / 2));

  if (m_radius == m_param->GetMinRadius()) {
    m_ready_to_divide = true;  
  }

}

Cell* Cell::Divide() {

  double c1_x, c1_y, c2_x, c2_y;

  c1_x = m_coordinates.first - m_radius * cos(m_axis.second);
  c1_y = m_coordinates.second - m_radius * sin(m_axis.second);
  c2_x = m_coordinates.first + m_radius * cos(m_axis.second);
  c2_y = m_coordinates.second + m_radius * sin(m_axis.second);

  m_coordinates = std::make_pair(c1_x, c1_y);
  m_axis.first = 0;
  m_axis.second = 0;
  m_radius = m_param->GetMinRadius();
  m_ready_to_divide = false;
	m_in_mitosis = false;
  return new Cell(std::make_pair(c2_x, c2_y), m_param, m_growth_rate);

}

bool Cell::ReadyToDivide() {

  return m_ready_to_divide;

}

std::pair<double, double> Cell::GetCoord() {

  return m_coordinates;

}

double Cell::GetRadius() {

	return m_radius;

}

double Cell::GetAxisLength() {

  return m_axis.first;

}

double Cell::GetAxisAngle() {

  return m_axis.second;

}

void Cell::SetGrowth(double rate) {

  m_growth_rate = rate;

}

double Cell::GetGrowth() {

  return m_growth_rate;

}

double Cell::CellDistance(Cell &other) {

  std::vector<double> centers;
  
  double a_x = GetCoord().first;
  double a_y = GetCoord().second;
  double a_rad = GetRadius();
  double a_len = GetAxisLength();
  double a_ang = GetAxisAngle();

  double b_x = other.GetCoord().first;
  double b_y = other.GetCoord().second;
  double b_rad = other.GetRadius();
  double b_len = other.GetAxisLength();
  double b_ang = other.GetAxisAngle();

  if (a_len > 0) {  

    centers.push_back(a_x + (-0.5 * a_len + a_rad) * cos(a_ang));
    centers.push_back(a_y + (-0.5 * a_len + a_rad) * sin(a_ang));
    centers.push_back(a_x + (0.5 * a_len + a_rad) * cos(a_ang));
    centers.push_back(a_y + (0.5 * a_len + a_rad) * sin(a_ang));

  } else {

    centers.push_back(a_x);
    centers.push_back(a_y);

  }

  if (b_len > 0) {

    centers.push_back(b_x + (-0.5 * b_len + b_rad) * cos(b_ang));
    centers.push_back(b_y + (-0.5 * b_len + b_rad) * sin(b_ang));
    centers.push_back(b_x + (0.5 * b_len + b_rad) * cos(b_ang));
    centers.push_back(b_y + (0.5 * b_len + b_rad) * sin(b_ang));

  } else {

    centers.push_back(b_x);
    centers.push_back(b_y);

  }    
  
  double x_dist, y_dist;
  double dist, min_dist = std::numeric_limits<double>::max();

  for (unsigned int i = 0; i < centers.size() / 2 - 1; i++) {

    for (unsigned int j = i + 1; j < centers.size() / 2; j++) {

      x_dist = centers[2*j] - centers[2*i];
      y_dist = centers[2*j+1] - centers[2*i+1]; 
      dist = pow(pow(x_dist,2) + pow(y_dist,2), 0.5) - a_rad - b_rad;
      min_dist = std::min(min_dist, dist);

    }

  }

  return min_dist;

} 



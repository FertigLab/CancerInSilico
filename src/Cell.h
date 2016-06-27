#ifndef CELL_HPP
#define CELL_HPP

#include <utility>
#include <vector>

#include "Parameters.h" 

class Cell {

private:

  Parameters* m_param;

  std::pair<double, double> m_coordinates;
  double m_radius;
  bool m_ready_to_divide;
  bool m_in_mitosis;
  double m_growth_rate;
  double m_min_cell_dist;
  std::pair<double, double> m_axis; //(length, angle)

public:

  Cell(std::pair<double, double>, Parameters*);
  Cell(std::pair<double, double>, Parameters*, double);

  Cell& DoTrial();
  void Migration();
  void Growth();
  void Rotation();
  void Deformation();
  bool ReadyToDivide();

	std::pair<double, double> GetCoord();
	double GetRadius();
  double GetAxisLength();
  double GetAxisAngle();
  Cell* Divide();
  void SetGrowth(double);
  double GetGrowth();

  double CellDistance(Cell&);

  bool operator!=(const Cell& b) const {
    return m_coordinates != b.m_coordinates;
  }

  bool operator==(const Cell& b) const {
    return m_coordinates == b.m_coordinates;
  }


};
 
#endif

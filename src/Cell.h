#ifndef CELL_HPP
#define CELL_HPP

#include "Point.h"
#include "Parameters.h"

class Cell {

  private:

    Parameters *m_param;

    Point m_coordinates;
    double m_radius;
    bool m_ready_to_divide;
    bool m_in_mitosis;
    double m_growth_rate;
	double m_axis_len, m_axis_ang;

  public:

	Cell(const Cell&, double);
    Cell(Point, Parameters*);
	Cell(Point, Parameters*, double);

    bool DoTrial(); //return true if growth
    void Translation();
    void Growth();
    void Rotation();
    void Deformation();

    Point GetCoord() const;
    void SetCoord(Point);
    double GetRadius() const;
	void SetRadius(double);
	void SetAxisLength(double);
    double GetAxisLength() const;
    double GetAxisAngle() const;
    void SetGrowth(double);
    double GetGrowth() const;
    bool ReadyToDivide() const;
	void EnterRandomPointOfMitosis();

 	Cell Divide();

    double CellDistance(const Cell&) const;

    bool operator!=(const Cell& b) const {
        return m_coordinates != b.m_coordinates;
    }

    bool operator==(const Cell& b) const {
        return m_coordinates == b.m_coordinates;
    }

    void operator=(const Cell& other) {

		m_param = other.m_param;
		m_coordinates = other.m_coordinates;
		m_in_mitosis = other.m_in_mitosis;
		m_ready_to_divide = other.m_ready_to_divide;
		m_axis_len = other.m_axis_len;
		m_axis_ang = other.m_axis_ang;
		m_radius = other.m_radius;
		m_growth_rate = other.m_growth_rate;

    }

};

#endif

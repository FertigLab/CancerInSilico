#include <cmath>
#include <Rcpp.h>

#include "Cell.h"

//used only for the initial population of cells
Cell::Cell(Point coord, Parameters* par) {

    m_param = par;
    m_coordinates = coord;
    m_in_mitosis = false;
    m_ready_to_divide = false;
    m_axis_ang = R::runif(0,2 * M_PI);
    m_radius = 1;
    m_axis_len = 2 * m_radius;
    m_growth_rate = 0;

}

//used only for daughter cells
Cell::Cell(Point coord, Parameters* par, double gr_rate) {

    m_param = par;
    m_coordinates = coord;
    m_in_mitosis = false;
    m_ready_to_divide = false;
    m_axis_len = 2;
    m_axis_ang = R::runif(0,2 * M_PI);
    m_radius = 1;
    m_growth_rate = gr_rate;

}

//should only be called for daughter cells (defn probably unneccesary)
Cell::Cell(const Cell& other, double gr_rate) {

    m_param = other.m_param;
    m_coordinates = other.m_coordinates;
    m_in_mitosis = other.m_in_mitosis;
    m_ready_to_divide = other.m_ready_to_divide;
    m_axis_len = other.m_axis_len;
    m_axis_ang = other.m_axis_ang;
    m_radius = other.m_radius;
    m_growth_rate = gr_rate;

}

Cell Cell::Divide() {

    double x = m_coordinates.x - cos(m_axis_ang);
    double y = m_coordinates.y - sin(m_axis_ang);

    m_coordinates.x += cos(m_axis_ang);
    m_coordinates.y += sin(m_axis_ang);
    m_axis_len = 2;
    m_axis_ang = 0;
    m_radius = 1;
    m_ready_to_divide = false;
    m_in_mitosis = false;

    return Cell(Point(x,y), m_param, m_growth_rate);

}

bool Cell::DoTrial() {

    double unif = R::runif(0, 1);
    double nG = m_param->GetNG();
    bool growth = false;

    if (!m_in_mitosis) { //Interphase

        if (unif <= (1.0 / (nG + 1.0))) { growth = true; Growth();}
        else { Translation();}

    } else { //Mitosis

        if (unif <= (1.0 / (nG + 1.0))) { growth = true; Deformation();}
        else if ((nG + 1.0) * unif <= 1.0 + nG / 2.0) { Rotation();}
        else if (!m_ready_to_divide) { Translation();}

    }

    return growth;

}

void Cell::Translation() {

    double length = m_param->GetMaxTranslation() * pow(R::runif(0, 1), 0.5);
    double direction = R::runif(0, 2 * M_PI);
    m_coordinates.x += length * cos(direction);
    m_coordinates.y += length * sin(direction);

}

void Cell::Growth() {

    double growth = R::runif(0,m_growth_rate);
    m_radius = std::min(pow(2,0.5), m_radius + growth);
    m_axis_len = 2 * m_radius;

    if (m_radius == pow(2,0.5)) {

        m_axis_ang = R::runif(0,2 * M_PI);
        m_in_mitosis = true;

    }

}

void Cell::Rotation() {

    double rotate = R::runif(-m_param->GetMaxRotate(), m_param->GetMaxRotate());
    m_axis_ang += rotate;

}


void Cell::Deformation() {

    double deform = R::runif(0, m_param->GetMaxDeform());
    m_axis_len = std::min(4.0, m_axis_len + deform);
    m_radius = m_param->GetRadius(m_axis_len);

    if (m_axis_len == 4.0) {

        m_ready_to_divide = true;

    }

}

bool Cell::ReadyToDivide() const {

    return m_ready_to_divide;

}

Point Cell::GetCoord() const {

    return m_coordinates;

}

void Cell::SetCoord(Point pt) {
    
    m_coordinates = pt;

}

double Cell::GetRadius() const {

    return m_radius;

}

void Cell::SetRadius(double rad) {

    m_radius = rad;

}

void Cell::SetAxisLength(double len) {

    m_axis_len = len;

}

double Cell::GetAxisLength() const {

    return m_axis_len;

}

double Cell::GetAxisAngle() const {

    return m_axis_ang;

}

void Cell::SetGrowth(double rate) {

    m_growth_rate = rate;

}

double Cell::GetGrowth() const {

    return m_growth_rate;

}

int Cell::GetCellType() const {

    return m_cell_type;

}

void Cell::SetCellType(int newCellType) {
    m_cell_type = newCellType;
}

void Cell::EnterRandomPointOfMitosis() {

    m_in_mitosis = true;
    m_axis_len = R::runif(2 * pow(2,0.5),4);
    m_radius = m_param->GetRadius(m_axis_len);
    
}

double Cell::GetArea() const {

    if (!m_in_mitosis) {
    
        return M_PI * pow(m_radius, 2);
  
    } else {

        return M_PI * pow(m_param->GetMaxRadius(), 2);

    }

}

double Cell::CellDistance(const Cell& other) const {

    std::vector<Point> a_centers;
    std::vector<Point> b_centers;

    double a_x = GetCoord().x;
    double a_y = GetCoord().y;
    double a_rad = GetRadius();
    double a_len = GetAxisLength();
    double a_ang = GetAxisAngle();

    double b_x = other.GetCoord().x;
    double b_y = other.GetCoord().y;
    double b_rad = other.GetRadius();
    double b_len = other.GetAxisLength();
    double b_ang = other.GetAxisAngle();

    a_centers.push_back(Point(
        a_x + (0.5 * a_len - a_rad) * cos(a_ang),
        a_y + (0.5 * a_len - a_rad) * sin(a_ang)));

    a_centers.push_back(Point(
        a_x - (0.5 * a_len - a_rad) * cos(a_ang),
        a_y - (0.5 * a_len - a_rad) * sin(a_ang)));

    b_centers.push_back(Point(
        b_x + (0.5 * b_len - b_rad) * cos(b_ang),
        b_y + (0.5 * b_len - b_rad) * sin(b_ang)));

    b_centers.push_back(Point(
        b_x - (0.5 * b_len - b_rad) * cos(b_ang),
        b_y - (0.5 * b_len - b_rad) * sin(b_ang)));

    double min_dist = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < a_centers.size(); ++i) {

        for (unsigned int j = 0; j < b_centers.size(); ++j) {

            min_dist = std::min(min_dist, a_centers[i].dist(b_centers[j]));

        }

    }

    return min_dist - a_rad - b_rad;

}



#include <cmath>
#include <Rcpp.h>

#include "Cell.h"

//used only for the initial population of cells
Cell::Cell(Point coord, Parameters* par) {

    mParams = par;
    mCoordinates = coord;
    mInMitosis = false;
    mReadyToDivide = false;
    mAxisAng = R::runif(0,2 * M_PI);
    mRadius = 1;
    mAxisLen = 2 * mRadius;
    mGrowthRate = 0;

}

//used only for daughter cells
Cell::Cell(Point coord, Parameters* par, double gr_rate) {

    mParams = par;
    mCoordinates = coord;
    mInMitosis = false;
    mReadyToDivide = false;
    mAxisLen = 2;
    mAxisAng = R::runif(0,2 * M_PI);
    mRadius = 1;
    mGrowthRate = gr_rate;

}

//should only be called for daughter cells (defn probably unneccesary)
Cell::Cell(const Cell& other, double gr_rate) {

    mParams = other.mParams;
    mCoordinates = other.mCoordinates;
    mInMitosis = other.mInMitosis;
    mReadyToDivide = other.mReadyToDivide;
    mAxisLen = other.mAxisLen;
    mAxisAng = other.mAxisAng;
    mRadius = other.mRadius;
    mGrowthRate = gr_rate;

}

Cell Cell::Divide() {

    double x = mCoordinates.x - cos(mAxisAng);
    double y = mCoordinates.y - sin(mAxisAng);

    mCoordinates.x += cos(mAxisAng);
    mCoordinates.y += sin(mAxisAng);
    mAxisLen = 2;
    mAxisAng = 0;
    mRadius = 1;
    mReadyToDivide = false;
    mInMitosis = false;

    return Cell(Point(x,y), mParams, mGrowthRate);

}

bool Cell::DoTrial() {

    double unif = R::runif(0, 1);
    double nG = mParams->nG();
    bool growth = false;

    if (!mInMitosis) { //Interphase

        if (unif <= (1.0 / (nG + 1.0))) { growth = true; Growth();}
        else { Translation();}

    } else { //Mitosis

        if (unif <= (1.0 / (nG + 1.0))) { growth = true; Deformation();}
        else if ((nG + 1.0) * unif <= 1.0 + nG / 2.0) { Rotation();}
        else if (!mReadyToDivide) { Translation();}

    }

    return growth;

}

void Cell::Translation() {

    double length = mParams->maxTranslation() * pow(R::runif(0, 1), 0.5);
    double direction = R::runif(0, 2 * M_PI);
    mCoordinates.x += length * cos(direction);
    mCoordinates.y += length * sin(direction);

}

void Cell::Growth() {

    double growth = R::runif(0,mGrowthRate);
    mRadius = std::min(pow(2,0.5), mRadius + growth);
    mAxisLen = 2 * mRadius;

    if (mRadius == pow(2,0.5)) {

        mAxisAng = R::runif(0,2 * M_PI);
        mInMitosis = true;

    }

}

void Cell::Rotation() {

    double rotate = R::runif(-mParams->maxRotate(), mParams->maxRotate());
    mAxisAng += rotate;

}


void Cell::Deformation() {

    double deform = R::runif(0, mParams->maxDeform());
    mAxisLen = std::min(4.0, mAxisLen + deform);
    mRadius = mParams->GetRadius(mAxisLen);

    if (mAxisLen == 4.0) {

        mReadyToDivide = true;

    }

}

bool Cell::ReadyToDivide() const {

    return mReadyToDivide;

}

Point Cell::GetCoord() const {

    return mCoordinates;

}

void Cell::SetCoord(Point pt) {
    
    mCoordinates = pt;

}

double Cell::GetRadius() const {

    return mRadius;

}

void Cell::SetRadius(double rad) {

    mRadius = rad;

}

void Cell::SetAxisLength(double len) {

    mAxisLen = len;

}

double Cell::GetAxisLength() const {

    return mAxisLen;

}

double Cell::GetAxisAngle() const {

    return mAxisAng;

}

void Cell::SetGrowth(double rate) {

    mGrowthRate = rate;

}

double Cell::GetGrowth() const {

    return mGrowthRate;

}

void Cell::EnterRandomPointOfMitosis() {

    mInMitosis = true;
    mAxisLen = R::runif(2 * pow(2,0.5),4);
    mRadius = mParams->GetRadius(mAxisLen);
    
}

double Cell::GetArea() const {

    if (!mInMitosis) {
    
        return M_PI * pow(mRadius, 2);
  
    } else {

        return M_PI * pow(mParams->maxRadius(), 2);

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



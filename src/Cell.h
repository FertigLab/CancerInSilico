#ifndef CELL_HPP
#define CELL_HPP

#include "Point.h"
#include "Parameters.h"

class Cell {

  private:

    /* copy of parameters object */
    Parameters* mParams;

    Point mCoordinates;
    double mRadius;
    bool mReadyToDivide;
    bool mInMitosis;
    double mGrowthRate;
	double mAxisLen, mAxisAng;

    /* type of cell */
    char mType;

  public:

    Cell(Point, Parameters*);
	Cell(Point, const Cell&);

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
    double GetArea() const;
    char GetType() const;

    bool ReadyToDivide() const;
	void EnterRandomPointOfMitosis();

 	Cell Divide();

    double CellDistance(const Cell&) const;

    bool operator!=(const Cell& b) const {
        return mCoordinates != b.mCoordinates;
    }

    bool operator==(const Cell& b) const {
        return mCoordinates == b.mCoordinates;
    }

};

#endif

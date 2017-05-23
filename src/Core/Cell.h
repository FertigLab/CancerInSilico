#ifndef CIS_CELL_H
#define CIS_CELL_H

// top-level abstract class for a cell

#include <Rcpp.h>
#include <vector>
#include <stdint.h>

#include "Drug.h"
#include "CellType.h"
#include "CellPhase.h"

class Cell
{
protected:

    // info about the cell's type
    CellType mType;

    // length of the cell cycle in hours
    double mCycleLength;

    // phase in the cell cycle
    CellPhase mPhase;

    // whether or not a drug has been applied
    uint64_t mDrugApplied;

    // true if mitosis is complete
    bool mReadyToDivide;

public:

    // constructors
    Cell(CellType);

    // getters
    CellPhase phase() const {return mPhase;}
    CellType type() const {return mType;}
    double cycleLength() const {return mCycleLength;}
    bool drugApplied(unsigned i) {return (mDrugApplied >> i) & 1;}
    bool readyToDivide() {return mReadyToDivide;}    

    // setters
    void setPhase(CellPhase phase) {mPhase = phase;}
    void setCycleLength(double len) {mCycleLength = len;}
    void setCellType(CellType type) {mType = type;}
    void setReadyToDivide(bool b) {mReadyToDivide = b;}

    // go to random point in the cell cycle
    virtual void gotoRandomCyclePoint() = 0;

    // apply drug to cell
    virtual double applyDrug(const Drug&);
};

#endif

#ifndef CIS_CELL_HPP
#define CIS_CELL_HPP

#include <vector.h>
#include <Rcpp.h>

#include "Point.h"
#include "Parameters.h"

class Cell
{

private:

    /* R object with all info about the cell's type */
    Rcpp::S4* mCellType;    

    /* must be R object so that R functions can manipulate it */
    Rcpp::List mProperties;

    /* phase in the cell cycle */
    enum CellPhase {I, M, G0, G1, S, G2};
    CellPhase mPhase;

public:

    Cell(Rcpp::S4*, Rcpp::List, CellPhase);

    /* NOT COPY CONSTRUCTOR - used for daughter cell, pass pointer to parent */
    Cell(const Cell*);
    virtual void divide(Cell*, const Cell*) = 0;

    /* getters */
    CellPhase phase() const {return mPhase;}
    Rcpp::S4* cellType() const {return mCellType;}
    virtual double area() const = 0;
    
    /* setters */   
    void setPhase(CellPhase phase) {mPhase = phase;}
    void setCellTpye(Rcpp::S4* type) {mCellType = type;}

    /* operators */
    virtual bool operator!=(const Cell&) const = 0;
    virtual bool operator==(const Cell&) const = 0;
    virtual double distance(const Cell&) const = 0;

    /* go to random point in the cell cycle */
    virtual void gotoRandomCyclePoint() = 0;
};

#endif

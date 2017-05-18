#ifndef CIS_CELL_TYPE_H
#define CIS_CELL_TYPE_H

// contains all info about a specific cell type

#include <Rcpp.h>

class CellType
{
private:

    // name of this cell type
    std::string mName;

    // id of the cell
    unsigned mID;

    // relative size of cell
    double mSize;

    // minimum possible cycle length    
    double mMinCycle;

    // needed to get cycle length
    Rcpp::S4 mCellTypeClass;

public:

    // constructor
    CellType(unsigned, const Rcpp::S4&);

    // getters
    std::string name() const  {return mName;}
    unsigned id()      const  {return mID;}
    double size()      const  {return mSize;}
    double minCycle()  const  {return mMinCycle;}    

    // get random cycle length based on this type's distribution
    double cycleLength() const;
};

#endif

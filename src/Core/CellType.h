#ifndef CIS_CELL_TYPE_H
#define CIS_CELL_TYPE_H

#include <Rcpp.h>

class CellType
{
private:

    // name of this cell type
    std::string mName;

    // numerical id of this type
    // TODO: id's are unneccesary
    unsigned mID;

    // relative size of cell
    double mSize;

    // minimum possible cycle length    
    double mMinCycle;

    // needed to get cycle length
    Rcpp::S4 mCellTypeClass;

public:

    CellType(unsigned, const Rcpp::S4&);

    std::string name() const  {return mName;}
    unsigned id() const       {return mID;}
    double size() const       {return mSize;}
    double minCycle() const   {return mMinCycle;}    

    double cycleLength() const;
};

#endif

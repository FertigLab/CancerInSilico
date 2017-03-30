#ifndef CIS_CELL_TYPE_H
#define CIS_CELL_TYPE_H

class CellType
{
private:

    // numerical id of this type
    unsigned mID;

    // relative size of cell
    double mSize;

    // distribution of cycle lengths
    double mCycleLengths[1000];

    // whether or not daughter cells inherit cycle length
    bool mInheritCycle;

public:

    CellType(unsigned, Rcpp::S4);

    unsigned id() const {return mID;}
    double size() const {return mSize;}
    bool inheritCycle() const {return mInheritCycle;}
    
    double cycleLength() const;
};

#endif

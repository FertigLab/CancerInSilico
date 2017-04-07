#ifndef CIS_CELL_TYPE_H
#define CIS_CELL_TYPE_H

#define DIST_SIZE 1000

class CellType
{
private:

    // name of this cell type
    std::string mName;

    // numerical id of this type
    unsigned mID;

    // relative size of cell
    double mSize;

    // distribution of cycle length
    double mCycleLength[DIST_SIZE];

public:

    CellType(unsigned, const Rcpp::S4&);

    unsigned id() const {return mID;}
    double size() const {return mSize;}
    std::string name() const {return mName;}
    
    double cycleLength() const;
    double minCycleLength() const;
};

#endif

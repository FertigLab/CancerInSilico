#ifndef CIS_PARAMETERS_H
#define CIS_PARAMETERS_H

#include <vector>
#include <Rcpp.h>

#include "Drug.h"
#include "CellType.h"

class Parameters
{
protected:

    // list of parameters from R 
    Rcpp::List mParams;

    // drugs used in the simulation
    std::vector<Drug> mDrugs;

    // possible cell types
    std::vector<CellType> mCellTypes;

public:

    // constructor
    Parameters(Rcpp::List);

    // return all parameters
    Rcpp::List getRParameters() {return mParams;}

    // general model parameters
    double initialNum()         {return mParams["initialNum"];}
    double runTime()            {return mParams["runTime"];}
    double density()            {return mParams["density"];}
    double boundary()           {return mParams["boundary"];}
    double randSeed()           {return mParams["randSeed"];}
    double syncCycles()         {return mParams["syncCycles"];}
    double outputIncrement()    {return mParams["outputIncrement"];}
    double recordIncrement()    {return mParams["recordIncrement"];}
    double timeIncrement()      {return mParams["timeIncrement"];}

    // set the boundary to a numeric value
    void setBoundary(double b) {mParams["boundary"] = b;}

    // get a random cell type from initial distribution
    CellType* randomCellType();    
};

#endif

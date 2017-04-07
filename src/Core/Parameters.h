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
    double randSeed()           {return mParams["randSeed"];}
    bool syncCycles()           {return mParams["syncCycles"];}
    double outputIncrement()    {return mParams["outputIncrement"];}
    double recordIncrement()    {return mParams["recordIncrement"];}
    double timeIncrement()      {return mParams["timeIncrement"];}
    bool boundary()    {return Rcpp::as<double>(mParams["boundary"]) >= 0.0;}

    // iterator to drug list
    std::vector<Drug>::iterator drugsBegin() {return mDrugs.begin();}
    std::vector<Drug>::iterator drugsEnd() {return mDrugs.end();}

    // set the boundary to a numeric value
    void setBoundary(double b) {mParams["boundary"] = b;}

    // get a random cell type from initial distribution
    const CellType& randomCellType();
};

#endif

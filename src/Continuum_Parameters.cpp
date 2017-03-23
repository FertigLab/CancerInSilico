#if 0

#ifndef CIS_PARAMETERS_H
#define CIS_PARAMETERS_H

#include <vector>
#include <Rcpp.h>

#include "RadiusSolver.h"

class Parameters
{
private:

    /* list of parameters from R */
    Rcpp::List mParams;

    /* radius solver for cell geometry calculations, declared static here 
       to prevent initialization process from running more than once */
    static RadiusSolver mSolver;

    /* process parameters */
    void StoreTimeIncrement();   
    void StoreDrasdoParameters();

public:

    /* constructor */
    Parameters(Rcpp::List);

    /* return all parameters */
    Rcpp::List GetRparameters() { return mParams;}

    /* general model parameters */
    double initialNum()         {return mParams["initialNum"];}
    double runTime()            {return mParams["runTime"];}
    double density()            {return mParams["density"];}
    double boundary()           {return mParams["boundary"];}
    double randSeed()           {return mParams["randSeed"];}
    double syncCycles()         {return mParams["syncCycles"];}
    double outputIncrement()    {return mParams["outputIncrement"];}
    double recordIncrement()    {return mParams["recordIncrement"];}
    double timeIncrement()      {return mParams["timeIncrement"];}

    /* Drasdo specific parameters */
    double nG() { return mParams["nG"];}
    double epsilon() { return mParams["epsilon"];}
    double delta() { return mParams["delta"];}
    double maxDeform() { return mParams["maxDeform"];}
    double maxTranslation() { return mParams["maxTranslation"];}
    double maxRotate() { return mParams["maxRotate"];}
   
    /* set the boundary to a numeric value */
    void setBoundary(double b) { mParams["boundary"] = b;}

    /* get a random cell type from initial distribution */
    Rcpp::S4* randomCellType();    

    /* return radius given axis length, perserves area of dumbell */
    double getRadius(double a) { return mSolver.getRadius(a);}
};

#endif

#endif

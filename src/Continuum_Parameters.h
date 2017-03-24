#ifndef CIS_CONTINUUM_PARAMETERS_H
#define CIS_CONTINUUM_PARAMETERS_H

#include <vector>
#include <Rcpp.h>

#include "Parameters.h"
#include "RadiusSolver.h"

class ContinuumParameters : public Parameters
{
private:

    // radius solver for cell geometry calculations, declared static here 
    // to prevent initialization process from running more than once
    static RadiusSolver mSolver;

    // calculate time increment based on provided parameters
    void CalculateTimeIncrement();

public:

    /* constructor */
    ContinuumParameters(Rcpp::List p) : Parameters(p) 
        {CalculateTimeIncrement();}

    // get parameters
    double nG()             {return mParams["nG"];}
    double epsilon()        {return mParams["epsilon"];}
    double delta()          {return mParams["delta"];}
    double maxDeform()      {return mParams["maxDeform"];}
    double maxTranslation() {return mParams["maxTranslation"];}
    double maxRotate()      {return mParams["maxRotate"];}
   
    /* return radius given axis length, perserves area of dumbell */
    double getRadius(double a) {return mSolver.radius(a);}
};

#endif

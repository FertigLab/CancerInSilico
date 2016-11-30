// [[Rcpp::depends(BH)]]

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <testthat.h>
#include <Rcpp.h>

class Parameters {

  private:

    /* list of parameters from R */
    Rcpp::List mParams;

    /* max radius of cell in Drasdo model */
    mMaxRadius;

    /* lookup tables for radius-axis values */
	std::vector<double> mSlowSolver;
	std::vector<double> mFastSolver;

    /* distribution of growth rates */
	std::vector<double> mGrowthDist;

    /* function that applies effect of drug to a growth rate */
    Rcpp::Function mDrugEffect;

    /* initialize lookup tables for radius-axis values */
    void InitializeRadiusSolver();
	void InitSlowSolver();
	void InitFastSolver();

    /* get theta (radius) given axis length, use slow lookup */
	double GetThetaSlow(double); 

    /* hash axis length for fast lookup table */
	int HashAxisLength(double);


public:

    Parameters(double, Rcpp::List);

    /* general model parameters */
    double initialNum() { return mParams["initialNum"];}
    double runTime() { return mParams["runTime"];}
    double density() { return mParams["density"];}
    double inheritGrowth() { return mParams["inheritGrowth"];}
    double drugTime() { return mParams["drugTime"];}
    double boundary() { return mParams["boundary"];}
    double randSeed() { return mParams["randSeed"];}
    double syncCycles() { return mParams["syncCycles"];}
    double outputIncrement() { return mParams["outputIncrement"];}
    double recordIncrement() { return mParams["recordIncrement"];}
    double timeIncrement() { return mParams["timeIncrement"];}

    /* Drasdo specific parameters */
    double nG() { return mParams["nG"];}
    double epsilon() { return mParams["epsilon"];}
    double delta() { return mParams["delta"];}
    double maxDeform() { return mParams["maxGrowth"];}
    double maxTranslation() { return mParams["maxGrowth"];}
    double maxRotate() { return mParams["maxGrowth"];}
   
    double GetDrugEffect(double gr) { return mDrugEffect(gr);}
	double GetRandomGrowthRate();

    double GetRadius(double);

};

#endif

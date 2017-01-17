// [[Rcpp::depends(BH)]]

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <Rcpp.h>

struct GreaterThan {

    bool operator()(double a, double b) const {
        
        return a > b;

    }

};

class Parameters {

  private:

    /* list of parameters from R */
    Rcpp::List mParams;

    /* max radius of cell in Drasdo model */
    double mMaxRadius;

    /* lookup tables for radius-axis values */
	std::vector<double> mSlowSolver;
	std::vector<double> mFastSolver;

    /* initialize lookup tables for radius-axis values */
    void InitializeRadiusSolver();
	void InitSlowSolver();
	void InitFastSolver();

    /* hash axis length for fast lookup table */
	int HashAxisLength(double);

    /* process parameters */
    void StoreTimeIncrement();   
    void StoreUpdateParameters();

public:

    Parameters(double, Rcpp::List);

    Rcpp::List GetRparameters() { return mParams;}

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
    double maxDeform() { return mParams["maxDeform"];}
    double maxTranslation() { return mParams["maxTranslation"];}
    double maxRotate() { return mParams["maxRotate"];}
   
    double maxRadius() { return mMaxRadius;}

    void setBoundary(double b) { mParams["boundary"] = b;}

    double GetDrugEffect(double);
	double GetRandomGrowthRate(char);
    char GetRandomCellType();

    double GetRadius(double);

    /* get theta (radius) given axis length, use slow lookup */
	double GetThetaSlow(double); 

};

#endif

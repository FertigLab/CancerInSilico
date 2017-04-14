#include <Rcpp.h>

#include "../TestHeader.h"
#include "../../Core/Random.h"
#include "../../Core/Drug.h"
#include "../../Core/Parameters.h"
#include "../../OffLatticeModel/OffLatticeParameters.h"
#include "../../CellModels/DrasdoHohmeModel.h"

CATCH_TEST_CASE("Test Parameters")
{
    Random::setSeed(12344445);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model = pkgEnv.find("drugsInSystem");

    Parameters params (&model);

    CATCH_REQUIRE(params.initialNum() == 100);
    CATCH_REQUIRE(params.runTime() == 4);
    CATCH_REQUIRE(params.density() == 0.4);
    CATCH_REQUIRE(params.randSeed() == 0);
    CATCH_REQUIRE(params.outputIncrement() == 4.0);
    CATCH_REQUIRE(params.recordIncrement() == 0.1);
    CATCH_REQUIRE(params.timeIncrement() == Approx(0.0007).epsilon(0.01));
    CATCH_REQUIRE(params.boundary() == 1);
    CATCH_REQUIRE(params.syncCycles() == false); 

    params.setBoundary(100.0);
    CATCH_REQUIRE(Rcpp::as<double>(model.slot("boundary")) == 100.0);

    int sum = 0;
    for (unsigned i = 0; i < 1000; ++i)
    {
        sum += params.randomCellType().id();
    }
    CATCH_REQUIRE(sum == 500);

    CellType c0 = params.randomCellType();
    CellType c1 = params.randomCellType();

    CATCH_REQUIRE(c0.id() == 0);
    CATCH_REQUIRE(c1.id() == 1);

    CATCH_REQUIRE(c0.name() == "DEFAULT");
    CATCH_REQUIRE(c0.size() == 1.0);
    CATCH_REQUIRE(c0.minCycle() == 48);

    CATCH_REQUIRE(c1.name() == "DOUBLE_SIZE");
    CATCH_REQUIRE(c1.size() == 2.0);
    CATCH_REQUIRE(c1.minCycle() == 48);

    DrugIterator it = params.drugsBegin();
    unsigned v = 0;
    for (; it != params.drugsEnd(); ++it)
    {
        CATCH_REQUIRE((*it).id() == v++);
    }

    Drug d0 = *(params.drugsBegin() + 0);
    Drug d1 = *(params.drugsBegin() + 1);
    Drug d2 = *(params.drugsBegin() + 2);
    Drug d3 = *(params.drugsBegin() + 3);
    
    CATCH_REQUIRE(d0.id() == 0);
    CATCH_REQUIRE(d0.timeAdded() == 0);
    CATCH_REQUIRE(d1.id() == 1);
    CATCH_REQUIRE(d1.timeAdded() == 0);
    CATCH_REQUIRE(d2.id() == 2);
    CATCH_REQUIRE(d2.timeAdded() == 0);
    CATCH_REQUIRE(d3.id() == 3);
    CATCH_REQUIRE(d3.timeAdded() == 6);
}

CATCH_TEST_CASE("Test OffLatticeParameters")
{
    Random::setSeed(12344445);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model = pkgEnv.find("drugsInSystem");

    OffLatticeParameters params (&model);

    CATCH_REQUIRE(params.initialNum() == 100);
    CATCH_REQUIRE(params.runTime() == 4);
    CATCH_REQUIRE(params.density() == 0.4);
    CATCH_REQUIRE(params.randSeed() == 0);
    CATCH_REQUIRE(params.outputIncrement() == 4.0);
    CATCH_REQUIRE(params.recordIncrement() == 0.1);
    CATCH_REQUIRE(params.timeIncrement() == Approx(0.0007).epsilon(0.01));
    CATCH_REQUIRE(params.boundary() == 100);
    CATCH_REQUIRE(params.syncCycles() == false); 
    CATCH_REQUIRE(params.maxDeformation() == 0.1);
    CATCH_REQUIRE(params.maxTranslation() == 0.1);
    CATCH_REQUIRE(params.maxRotation() == Approx(0.309).epsilon(0.01));
    CATCH_REQUIRE(params.maxRadius() == 2.0);
    
    int sum = 0;
    for (unsigned i = 0; i < 1000; ++i)
    {
        sum += params.randomCellType().id();
    }
    CATCH_REQUIRE(sum == 500);

    CellType c0 = params.randomCellType();
    CellType c1 = params.randomCellType();

    CATCH_REQUIRE(c0.id() == 0);
    CATCH_REQUIRE(c1.id() == 1);

    CATCH_REQUIRE(c0.name() == "DEFAULT");
    CATCH_REQUIRE(c0.size() == 1.0);
    CATCH_REQUIRE(c0.minCycle() == 48);

    CATCH_REQUIRE(c1.name() == "DOUBLE_SIZE");
    CATCH_REQUIRE(c1.size() == 2.0);
    CATCH_REQUIRE(c1.minCycle() == 48);

    DrugIterator it = params.drugsBegin();
    unsigned v = 0;
    for (; it != params.drugsEnd(); ++it)
    {
        CATCH_REQUIRE((*it).id() == v++);
    }

    Drug d0 = *(params.drugsBegin() + 0);
    Drug d1 = *(params.drugsBegin() + 1);
    Drug d2 = *(params.drugsBegin() + 2);
    Drug d3 = *(params.drugsBegin() + 3);
    
    CATCH_REQUIRE(d0.id() == 0);
    CATCH_REQUIRE(d0.timeAdded() == 0);
    CATCH_REQUIRE(d1.id() == 1);
    CATCH_REQUIRE(d1.timeAdded() == 0);
    CATCH_REQUIRE(d2.id() == 2);
    CATCH_REQUIRE(d2.timeAdded() == 0);
    CATCH_REQUIRE(d3.id() == 3);
    CATCH_REQUIRE(d3.timeAdded() == 6);
}

CATCH_TEST_CASE("Test DrasdoHohmeParameters")
{
    Random::setSeed(12344445);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model = pkgEnv.find("drugsInSystem");

    DrasdoHohmeParameters params (&model);

    CATCH_REQUIRE(params.initialNum() == 100);
    CATCH_REQUIRE(params.runTime() == 4);
    CATCH_REQUIRE(params.density() == 0.4);
    CATCH_REQUIRE(params.randSeed() == 0);
    CATCH_REQUIRE(params.outputIncrement() == 4.0);
    CATCH_REQUIRE(params.recordIncrement() == 0.1);
    CATCH_REQUIRE(params.timeIncrement() == Approx(0.0007).epsilon(0.01));
    CATCH_REQUIRE(params.boundary() == 100);
    CATCH_REQUIRE(params.syncCycles() == false); 
    CATCH_REQUIRE(params.maxDeformation() == 0.1);
    CATCH_REQUIRE(params.maxTranslation() == 0.1);
    CATCH_REQUIRE(params.maxRotation() == Approx(0.309).epsilon(0.01));
    CATCH_REQUIRE(params.maxRadius() == 2.0);
    CATCH_REQUIRE(params.nG() == 28);
    CATCH_REQUIRE(params.epsilon() == 10.0);
    CATCH_REQUIRE(params.delta() == 0.2);
    
    int sum = 0;
    for (unsigned i = 0; i < 1000; ++i)
    {
        sum += params.randomCellType().id();
    }
    CATCH_REQUIRE(sum == 500);

    CellType c0 = params.randomCellType();
    CellType c1 = params.randomCellType();

    CATCH_REQUIRE(c0.id() == 0);
    CATCH_REQUIRE(c1.id() == 1);

    CATCH_REQUIRE(c0.name() == "DEFAULT");
    CATCH_REQUIRE(c0.size() == 1.0);
    CATCH_REQUIRE(c0.minCycle() == 48);

    CATCH_REQUIRE(c1.name() == "DOUBLE_SIZE");
    CATCH_REQUIRE(c1.size() == 2.0);
    CATCH_REQUIRE(c1.minCycle() == 48);

    DrugIterator it = params.drugsBegin();
    unsigned v = 0;
    for (; it != params.drugsEnd(); ++it)
    {
        CATCH_REQUIRE((*it).id() == v++);
    }

    Drug d0 = *(params.drugsBegin() + 0);
    Drug d1 = *(params.drugsBegin() + 1);
    Drug d2 = *(params.drugsBegin() + 2);
    Drug d3 = *(params.drugsBegin() + 3);
    
    CATCH_REQUIRE(d0.id() == 0);
    CATCH_REQUIRE(d0.timeAdded() == 0);
    CATCH_REQUIRE(d1.id() == 1);
    CATCH_REQUIRE(d1.timeAdded() == 0);
    CATCH_REQUIRE(d2.id() == 2);
    CATCH_REQUIRE(d2.timeAdded() == 0);
    CATCH_REQUIRE(d3.id() == 3);
    CATCH_REQUIRE(d3.timeAdded() == 6);
}

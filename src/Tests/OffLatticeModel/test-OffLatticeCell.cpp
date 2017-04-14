#include <Rcpp.h>

#include "../TestHeader.h"
#include "../../Core/Random.h"
#include "../../OffLatticeModel/OffLatticeCell.h"

CATCH_TEST_CASE("Test Cell")
{
    Random::setSeed(0);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model1 = pkgEnv.find("manyCellTypes");
    Rcpp::S4 model2 = pkgEnv.find("drugsInSystem");

    Rcpp::List types = model1.slot("cellTypes");
    Rcpp::List drugs = model2.slot("drugs");

    Drug drug0 = Drug(0, drugs(0)); // NO_EFFECT
    Drug drug1 = Drug(1, drugs(1)); // HALF_CYCLE_LENGTH
    Drug drug2 = Drug(2, drugs(2)); // HALF_DEFAULT_TYPE
    Drug drug3 = Drug(3, drugs(3)); // ADD_LATE

    CellType type0 (0, types(0)); 
    CellType type1 (1, types(1)); 

    OffLatticeCell cell0 (type0);
    OffLatticeCell cell1 (type0);
    OffLatticeCell cell2 (type0);

    CATCH_REQUIRE(cell0.phase() == INTERPHASE);
    CATCH_REQUIRE(cell1.phase() == INTERPHASE);
    CATCH_REQUIRE(cell2.phase() == INTERPHASE);

    CATCH_REQUIRE(cell0.type().id() == 0);
    CATCH_REQUIRE(cell0.type().size() == 1);
    CATCH_REQUIRE(cell0.type().minCycle() == 48);

    CATCH_REQUIRE(cell0.cycleLength() == 48);
    CATCH_REQUIRE(cell1.cycleLength() == 48);
    CATCH_REQUIRE(cell2.cycleLength() == 48);

    CATCH_REQUIRE(cell0.cycleLength() == 48);
    CATCH_REQUIRE(cell1.cycleLength() == 48);
    CATCH_REQUIRE(cell2.cycleLength() == 48);

    CATCH_REQUIRE(!cell0.drugApplied(0));
    CATCH_REQUIRE(!cell0.drugApplied(1));
    CATCH_REQUIRE(!cell0.drugApplied(2));
    CATCH_REQUIRE(!cell0.drugApplied(3));

    CATCH_REQUIRE(!cell0.readyToDivide());
    CATCH_REQUIRE(!cell1.readyToDivide());
    CATCH_REQUIRE(!cell2.readyToDivide());

    cell0.setPhase(MITOSIS);
    cell0.setCycleLength(4.0);
    cell0.setCellType(type1);
    cell0.setReadyToDivide(true);

    CATCH_REQUIRE(cell0.phase() == MITOSIS);    
    CATCH_REQUIRE(cell0.cycleLength() == 4.0);
    CATCH_REQUIRE(cell0.type().id() == 1);
    CATCH_REQUIRE(cell0.readyToDivide());

    cell0.applyDrug(drug0);
    cell0.applyDrug(drug1);
    cell0.applyDrug(drug3);

    cell1.applyDrug(drug0);
    cell1.applyDrug(drug1);
    cell1.applyDrug(drug2);
    cell1.applyDrug(drug3);

    CATCH_REQUIRE(cell0.drugApplied(0));
    CATCH_REQUIRE(cell0.drugApplied(1));
    CATCH_REQUIRE(!cell0.drugApplied(2));
    CATCH_REQUIRE(cell0.drugApplied(3));
    CATCH_REQUIRE(!cell0.drugApplied(4));

    CATCH_REQUIRE(cell0.cycleLength() == 2.0);
}

CATCH_TEST_CASE("Test OffLatticeCell")
{
    Random::setSeed(0);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model1 = pkgEnv.find("manyCellTypes");
    Rcpp::S4 model2 = pkgEnv.find("drugsInSystem");

    Rcpp::List types = model1.slot("cellTypes");
    Rcpp::List drugs = model2.slot("drugs");

    Drug drug0 = Drug(0, drugs(0)); // NO_EFFECT
    Drug drug1 = Drug(1, drugs(1)); // HALF_CYCLE_LENGTH
    Drug drug2 = Drug(2, drugs(2)); // HALF_DEFAULT_TYPE
    Drug drug3 = Drug(3, drugs(3)); // ADD_LATE

    CellType type0 (0, types(0)); 
    CellType type1 (1, types(1)); 

    OffLatticeCell cell0 (type0);
    OffLatticeCell cell1 (type1);

    // test inherited properties
    CATCH_REQUIRE(cell0.phase() == INTERPHASE);
    CATCH_REQUIRE(cell0.type().id() == 0);
    CATCH_REQUIRE(cell0.type().size() == 1);
    CATCH_REQUIRE(cell0.type().minCycle() == 48);
    CATCH_REQUIRE(cell0.cycleLength() == 48);
    CATCH_REQUIRE(cell0.cycleLength() == 48);
    CATCH_REQUIRE(!cell0.drugApplied(0));
    CATCH_REQUIRE(!cell0.drugApplied(1));
    CATCH_REQUIRE(!cell0.drugApplied(2));
    CATCH_REQUIRE(!cell0.drugApplied(3));
    CATCH_REQUIRE(!cell0.readyToDivide());

    // test properties
    CATCH_REQUIRE(cell0.coordinates() == Point<double>(0,0));
    CATCH_REQUIRE(cell0.radius() == 1);
    CATCH_REQUIRE(cell1.radius() == sqrt(2));
    CATCH_REQUIRE(cell0.axisLength() == 2);
    CATCH_REQUIRE(cell1.axisLength() == 2 * sqrt(2));
    CATCH_REQUIRE(cell0.axisAngle() == Approx(3.45).epsilon(0.01));
    CATCH_REQUIRE(cell0.area() == Approx(M_PI));
    CATCH_REQUIRE(cell1.area() == Approx(2 * M_PI));

    CATCH_REQUIRE(cell0.centers().first == cell0.centers().second);
    CATCH_REQUIRE(cell1.centers().first == cell1.centers().second);

    CATCH_REQUIRE(cell0 == cell1);
    double dist = -cell0.radius() - cell1.radius();
    CATCH_REQUIRE(cell0.distance(cell1) == Approx(dist));
    CATCH_REQUIRE(cell1.distance(cell0) == Approx(dist));

    cell0.setCoordinates(Point<double>(3.0, 5.0));
    cell1.setCoordinates(Point<double>(9.0, -3.0));

    dist = 10.0 - cell0.radius() - cell1.radius();
    CATCH_REQUIRE(cell0.distance(cell1) == Approx(dist));
    CATCH_REQUIRE(cell1.distance(cell0) == Approx(dist));
    CATCH_REQUIRE(cell0 != cell1);
    
    cell0.setAxisLength(3.0);
    CATCH_REQUIRE(cell0.radius() == Approx(1.284).epsilon(0.01));

    cell1.setRadius(2.0);
    CATCH_REQUIRE(cell1.axisLength() == 4.0);

    for (unsigned i = 0; i < 1000; ++i)
    {
        cell0.gotoRandomCyclePoint();
        if (cell0.phase() == INTERPHASE)
        {
            CATCH_REQUIRE(cell0.radius() * 2 == cell0.axisLength());
            CATCH_REQUIRE(cell0.radius() >= 1);
            CATCH_REQUIRE(cell0.radius() <= sqrt(2));
        }
        else
        {
            CATCH_REQUIRE(cell0.axisLength() >= sqrt(2));
            CATCH_REQUIRE(cell0.axisLength() <= 4);
        }
    }

    cell0.setAxisLength(4.0);
    OffLatticeCell cell0_D (cell0.type());
    cell0.divide(cell0_D);

    CATCH_REQUIRE(cell0.coordinates().x == Approx(3.95).epsilon(0.01));
    CATCH_REQUIRE(cell0.coordinates().y == Approx(5.30).epsilon(0.01));
    CATCH_REQUIRE(cell0_D.coordinates().x == Approx(2.05).epsilon(0.01));
    CATCH_REQUIRE(cell0_D.coordinates().y == Approx(4.70).epsilon(0.01));

    CATCH_REQUIRE(cell0.type().id() == cell0_D.type().id());
    CATCH_REQUIRE(cell0.radius() == 1);
    CATCH_REQUIRE(cell0_D.radius() == 1);
    CATCH_REQUIRE(cell0.axisLength() == 2);
    CATCH_REQUIRE(cell0_D.axisLength() == 2);

    CATCH_REQUIRE(cell0.axisAngle() == Approx(2.30).epsilon(0.01));
    CATCH_REQUIRE(cell0_D.axisAngle() == Approx(0.06).epsilon(0.01));
}



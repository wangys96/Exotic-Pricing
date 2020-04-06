/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif
#include <ql/termstructures/yield/all.hpp>
#include <ql/instruments/asianoption.hpp>
#include <ql/pricingengines/asian/all.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/utilities/dataformatters.hpp>

#include <iostream>
#include <iomanip>

using namespace QuantLib;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {
    Integer sessionId() { return 0; }
}
#endif

namespace {
    std::string averageTypeToString(Average::Type averageType) {
        if (averageType == Average::Geometric)
            return "Geometric Averaging";
        else if (averageType == Average::Arithmetic)
            return "Arithmetic Averaging";
        else
            return "UNKNOWN";
    }}

int main(int, char* []) {

    try {
        std::cout << std::endl;

        // set up dates
        Calendar calendar = TARGET();
        Date todaysDate(3, Jun, 2019);
        Date startDate(2, Jul, 2019);
        Date maturity(3, Sep, 2019);
        Settings::instance().evaluationDate() = todaysDate;

        // our options
        Average::Type avgType = Average::Geometric;
        Option::Type optType = Option::Call;
        Real underlying = 1;
        Real strike = 1;
        Rate riskFreeRate = 0.027238;
        Volatility volatility = 0.34054223610725;
       
        DayCounter dayCounter = Actual365Fixed();

        std::cout << "Average type = " << averageTypeToString(avgType) << std::endl;
        std::cout << "Option type = " << optType << std::endl;
        std::cout << "Today = " << todaysDate << std::endl;
        std::cout << "Fixing start = " << startDate << std::endl;
        std::cout << "Maturity = " << maturity << std::endl;
        std::cout << "Underlying = " << underlying << std::endl;
        std::cout << "Strike = " << strike << std::endl;;
        std::cout << std::endl;
        std::string method;
        std::cout << std::endl;

        // write column headings
        /*Size widths[] = { 35, 14, 14, 14 };
        std::cout << std::setw(widths[0]) << std::left << "Method"
            << std::setw(widths[1]) << std::left << "European"
            << std::setw(widths[2]) << std::left << "Bermudan"
            << std::setw(widths[3]) << std::left << "American"
            << std::endl;*/

        std::vector<Date> fixingDates;
        for (Integer i = 0; i <= 63; i++)
            fixingDates.push_back(startDate + i*Days);

        Handle<Quote> underlyingH(
            ext::shared_ptr<Quote>(new SimpleQuote(underlying)));

        // bootstrap the yield/dividend/vol curves    
        Handle<YieldTermStructure> flatTermStructure(
            ext::shared_ptr<YieldTermStructure>(
                new FlatForward(todaysDate, riskFreeRate, dayCounter)));

        Handle<BlackVolTermStructure> flatVolTS(
            ext::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(todaysDate, calendar, volatility,
                    dayCounter)));

        ext::shared_ptr<BlackScholesProcess> bsProcess(
            new BlackScholesProcess(underlyingH, flatTermStructure, flatVolTS));

        ext::shared_ptr<StrikedTypePayoff> payoff(
            new PlainVanillaPayoff(optType, strike));

        ext::shared_ptr<Exercise> europeanExercise(
            new EuropeanExercise(maturity));

        // options
        DiscreteAveragingAsianOption asianOption(
            avgType, 1.0, 0,fixingDates, payoff, europeanExercise);

        // Analytic for Asian
        method = "Analytic";
        asianOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
            new AnalyticDiscreteGeometricAveragePriceAsianEngine(bsProcess)));
        std::cout << "AP Analytic NPV = " << asianOption.NPV() << std::endl;
        asianOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
            new AnalyticDiscreteGeometricAverageStrikeAsianEngine(bsProcess)));
        std::cout << "AS Analytic NPV = " << asianOption.NPV() << std::endl;

        // Monte Carlo Method
        method = "QMC";
        Size mcSeed = 42;
        ext::shared_ptr<PricingEngine> mcengine1;
        mcengine1 = MakeMCDiscreteGeometricAPEngine<LowDiscrepancy>(bsProcess)
            .withSamples(250000)
            //.withAbsoluteTolerance(0.0001)
            .withBrownianBridge(true)
            .withSeed(mcSeed);
        asianOption.setPricingEngine(mcengine1);
        std::cout << "AP QMC NPV = " << asianOption.NPV() << std::endl;

        ext::shared_ptr<PricingEngine> mcengine2;
        mcengine2 = MakeMCDiscreteGeometricASEngine<LowDiscrepancy>(bsProcess)
            .withSamples(250000)
            //.withAbsoluteTolerance(0.0001)
            .withBrownianBridge(true)
            .withSeed(mcSeed);
        asianOption.setPricingEngine(mcengine2);
        std::cout << "AS QMC NPV = " << asianOption.NPV() << std::endl;


        // End test
        system("PAUSE");
        return 0;

    }
    catch (std::exception & e) {
        std::cerr << e.what() << std::endl;
        system("PAUSE");
        return 1;
    }
    catch (...) {
        std::cerr << "unknown error" << std::endl;
        system("PAUSE");
        return 1;
    }
}
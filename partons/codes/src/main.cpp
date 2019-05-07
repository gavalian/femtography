#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <ElementaryUtils/parameters/Parameter.h>
#include <ElementaryUtils/string_utils/Formatter.h>
#include <partons/beans/convol_coeff_function/DVCS/DVCSConvolCoeffFunctionKinematic.h>
#include <partons/beans/gpd/GPDKinematic.h>
#include <partons/beans/gpd/GPDType.h>
#include <partons/beans/observable/ObservableKinematic.h>
#include <partons/beans/parton_distribution/PartonDistribution.h>
#include <partons/beans/parton_distribution/QuarkDistribution.h>
#include <partons/beans/PerturbativeQCDOrderType.h>
#include <partons/beans/QuarkFlavor.h>
#include <partons/FundamentalPhysicalConstants.h>
#include <partons/modules/convol_coeff_function/DVCS/DVCSCFFStandard.h>
#include <partons/modules/gpd/GPDGK16Numerical.h>
#include <partons/modules/gpd/GPDMMS13.h>
#include <partons/modules/observable/DVCS/asymmetry/DVCSAulMinus.h>
#include <partons/modules/process/DVCS/DVCSProcessBMJ12.h>
#include <partons/modules/scales/ScalesQ2Multiplier.h>
#include <partons/modules/xi_converter/XiConverterXBToXi.h>
#include <partons/ModuleObjectFactory.h>
#include <partons/Partons.h>
#include <partons/services/ConvolCoeffFunctionService.h>
#include <partons/services/GPDService.h>
#include <partons/services/ObservableService.h>
#include <partons/ServiceObjectRegistry.h>
#include <QtCore/qcoreapplication.h>
#include <stddef.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

#include "../include/DVCSCrossSectionUUMinusBHSubProc.h"
#include "../include/DVCSCrossSectionUUMinusDVCSSubProc.h"
#include "../include/DVCSCrossSectionUUMinusINTSubProc.h"

// Typical value
// Number of steps
// Range

// x
const double c_X = 0.1;
const size_t c_nX = 100;
const std::pair<double, double> c_rangeX = std::make_pair(-0.95, 0.95);

// xi
const double c_Xi = 0.2;
const size_t c_nXi = 100;
const std::pair<double, double> c_rangeXi = std::make_pair(1.E-6, 0.95);

// t
const double c_TAbs = 0.3;
const size_t c_nTAbs = 100;
const std::pair<double, double> c_rangeTAbs = std::make_pair(0.1, 0.5);

// xB
const double c_Xb = 0.4;
const size_t c_nXb = 100;
const std::pair<double, double> c_rangeXb = std::make_pair(1.E-6, 0.95);

// Q2 (we assume Q2 = muF2 = muR2)
const double c_Q2 = 2.;
const size_t c_nQ2 = 100;
const std::pair<double, double> c_rangeQ2 = std::make_pair(1., 10.);

// E
const double c_E = 11.;
const size_t c_nE = 10;
const std::pair<double, double> c_rangeE = std::make_pair(1., 200.);

// phi
const double c_Phi = M_PI / 2.;
const size_t c_nPhi = 10;
const std::pair<double, double> c_rangePhi = std::make_pair(0., 2 * M_PI);

/*
 * Evaluate value on a loop
 */

double getValue(const size_t i, const size_t n,
        const std::pair<double, double>& range, const bool isLog) {

    if (isLog) {
        /*double step = (log10(range.second) - log10(range.first))
                      / double(n);
        double value = step*0.5 + i * step;
        return pow(10.0,value);*/
        return pow(10.,
                log10(range.first)
                        + (i+1) * (log10(range.second) - log10(range.first))
                                / double(n));
    } else {
      // Changed by Gagik.
        double step = (range.second - range.first) / double(n);
        return 0.5*step + i*step;
        //return range.first + i * (range.second - range.first) / double(n);
    }
}

bool isDVCSKinematicsValid(const double xB, const double t, const double Q2,
        const double E) {

    double m_epsilon = 2 * xB * PARTONS::Constant::PROTON_MASS / sqrt(Q2);
    double m_y = Q2 / (2 * xB * PARTONS::Constant::PROTON_MASS * E);
    double eps2 = m_epsilon * m_epsilon;
    double epsroot = sqrt(1. + eps2);
    double tfactor = -Q2 / (4 * xB * (1. - xB) + eps2);
    double m_tmin = tfactor * (2 * (1. - xB) * (1 - epsroot) + eps2);
    double m_tmax = tfactor * (2 * (1. - xB) * (1 + epsroot) + eps2);
    double m_xBmin = 2 * Q2 * E / PARTONS::Constant::PROTON_MASS
            / (4 * E * E - Q2);

    if (xB < m_xBmin || xB > 1.)
        return false;
    if (t > m_tmin || t < m_tmax)
        return false;
    if (Q2 < 0.)
        return false;
    if (E < 0.)
        return false;
    if (m_y < 0. || m_y > 1.)
        return false;

    return true;
}

/*
 * Main function.
 */
int main(int argc, char** argv) {

    // Init Qt4
    QCoreApplication a(argc, argv);
    PARTONS::Partons* pPartons = 0;

    // Init PARTONS application
    pPartons = PARTONS::Partons::getInstance();
    pPartons->init(argc, argv);

    try {

        // Mode
        if (argc != 3) {
            throw ElemUtils::CustomException("main", __func__,
                    ElemUtils::Formatter() << "Usage: " << argv[0]
                            << " mode output_file_path");
        }

        if (std::string(argv[1]) != "GPDGK16Numerical"
                && std::string(argv[1]) != "GPDMMS13"
                && std::string(argv[1]) != "CFF"
                && std::string(argv[1]) != "OBS_ALU"
                && std::string(argv[1]) != "OBS_CS") {
            throw ElemUtils::CustomException("main", __func__,
                    ElemUtils::Formatter() << "Unknown mode: " << argv[1]);
        }

        if (std::string(argv[1]) == "GPDGK16Numerical"
                || std::string(argv[1]) == "GPDMMS13") {
            PARTONS::Partons::getInstance()->getLoggerManager()->info("main",
                    __func__,
                    ElemUtils::Formatter() << "GPD H for up quarks and "
                            << argv[1]
                            << " model will be evaluated as a function of (x, xi, -t) for muF2 = muR2 = "
                            << c_Q2 << " GeV2");
        } else if (std::string(argv[1]) == "CFF") {
            PARTONS::Partons::getInstance()->getLoggerManager()->info("main",
                    __func__,
                    ElemUtils::Formatter()
                            << "DVCS LO imaginary CFF H will be evaluated as a function of (xi, -t, Q2) for GPD model GPDGK16Numerical");
        } else if (std::string(argv[1]) == "OBS_ALU") {
            PARTONS::Partons::getInstance()->getLoggerManager()->info("main",
                    __func__,
                    ElemUtils::Formatter()
                            << "DVCS ALU observable will be evaluated as a function of (xB, -t, phi) for GPD model GPDGK16Numerical and LO CFF and Q2 = "
                            << c_Q2 << " GeV2 and E = " << c_E << " GeV");
        } else if (std::string(argv[1]) == "OBS_CS") {
            PARTONS::Partons::getInstance()->getLoggerManager()->info("main",
                    __func__,
                    ElemUtils::Formatter()
                            << "DVCS sigmaUU observable will be evaluated as a function of (xB, phi, sub_process) for GPD model GPDGK16Numerical and LO CFF and t = "
                            << -1 * c_TAbs << " GeV2 and Q2 = " << c_Q2
                            << " GeV2 and E = " << c_E << " GeV");
        }

        // Retrieve GPD service
        PARTONS::GPDService* pGPDService =
                PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

        // Retrieve service
        PARTONS::ConvolCoeffFunctionService* pDVCSConvolCoeffFunctionService =
                PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getConvolCoeffFunctionService();

        // Retrieve Observable service
        PARTONS::ObservableService* pObservableService =
                PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getObservableService();

        // Create GPDModule
        PARTONS::GPDModule* pGPDModule;

        if (std::string(argv[1]) == "GPDMMS13") {
            pGPDModule =
                    PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                            PARTONS::GPDMMS13::classId);
        } else {
            pGPDModule =
                    PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                            PARTONS::GPDGK16Numerical::classId);
        }

        // Create CFF module
        PARTONS::DVCSConvolCoeffFunctionModule* pDVCSCFFModel =
                PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSConvolCoeffFunctionModule(
                        PARTONS::DVCSCFFStandard::classId);

        // Set its PerturbativeQCDOrder
        pDVCSCFFModel->configure(
                ElemUtils::Parameter(
                        PARTONS::PerturbativeQCDOrderType::PARAMETER_NAME_PERTURBATIVE_QCD_ORDER_TYPE,
                        PARTONS::PerturbativeQCDOrderType::LO));

        // Create XiConverterModule
        PARTONS::XiConverterModule* pXiConverterModule =
                PARTONS::Partons::getInstance()->getModuleObjectFactory()->newXiConverterModule(
                        PARTONS::XiConverterXBToXi::classId);

        // Create ScalesModule
        PARTONS::ScalesModule* pScalesModule =
                PARTONS::Partons::getInstance()->getModuleObjectFactory()->newScalesModule(
                        PARTONS::ScalesQ2Multiplier::classId);

        // Create ProcessModule
        PARTONS::DVCSProcessModule* pProcessModule =
                PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSProcessModule(
                        PARTONS::DVCSProcessBMJ12::classId);

        // Create Observable
        PARTONS::Observable* pObservable[3];

        for (size_t i = 0; i < 3; i++)
            pObservable[i] = 0;

        if (std::string(argv[1]) == "OBS_ALU") {

            pObservable[0] =
                    PARTONS::Partons::getInstance()->getModuleObjectFactory()->newObservable(
                            PARTONS::DVCSAulMinus::classId);

            pObservable[0]->setProcessModule(pProcessModule);
        }

        if (std::string(argv[1]) == "OBS_CS") {

            pObservable[0] =
                    PARTONS::Partons::getInstance()->getModuleObjectFactory()->newObservable(
                            PARTONS::DVCSCrossSectionUUMinusBHSubProc::classId);
            pObservable[1] =
                    PARTONS::Partons::getInstance()->getModuleObjectFactory()->newObservable(
                            PARTONS::DVCSCrossSectionUUMinusINTSubProc::classId);
            pObservable[2] =
                    PARTONS::Partons::getInstance()->getModuleObjectFactory()->newObservable(
                            PARTONS::DVCSCrossSectionUUMinusDVCSSubProc::classId);

            pObservable[0]->setProcessModule(pProcessModule);
            pObservable[1]->setProcessModule(pProcessModule);
            pObservable[2]->setProcessModule(pProcessModule);
        }

        // Link modules (set physics assumptions of your computation)

        pProcessModule->setScaleModule(pScalesModule);
        pProcessModule->setXiConverterModule(pXiConverterModule);
        pProcessModule->setConvolCoeffFunctionModule(pDVCSCFFModel);
        pDVCSCFFModel->setGPDModule(pGPDModule);

        // Open file
        std::ofstream outputFile;

        outputFile.open(argv[2]);
        outputFile << std::scientific;

        if (!outputFile.is_open()) {
            throw ElemUtils::CustomException("main", __func__,
                    ElemUtils::Formatter() << "Cannot open: " << argv[2]);
        }

        // *******************
        // Run GPD computation
        // *******************
        if (std::string(argv[1]) == "GPDGK16Numerical"
                || std::string(argv[1]) == "GPDMMS13") {

            //counter
            size_t counterI = 0;
            size_t counterN = (c_nX + 1) * (c_nXi + 1) * (c_nTAbs + 1);

            //loop over x
            for (size_t iX = 0; iX <= c_nX; iX++) {

                double thisX = getValue(iX, c_nX, c_rangeX, false);

                //loop over xi
                for (size_t iXi = 0; iXi <= c_nXi; iXi++) {

                    double thisXi = getValue(iXi, c_nXi, c_rangeXi, true);

                    //loop over t
                    for (size_t iTAbs = 0; iTAbs <= c_nTAbs; iTAbs++) {

                        double thisTAbs = getValue(iTAbs, c_nTAbs, c_rangeTAbs,
                                false);

                        // Run computation
                        PARTONS::GPDResult gpdResult =
                                pGPDService->computeGPDModel(
                                        PARTONS::GPDKinematic(thisX, thisXi,
                                                -1 * thisTAbs, c_Q2, c_Q2),
                                        pGPDModule);

                        // Write to file
                        outputFile << thisX << "\t" << thisXi << "\t"
                                << thisTAbs << "\t"
                                << gpdResult.getPartonDistribution(
                                        PARTONS::GPDType::H).getQuarkDistribution(
                                        PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                                << std::endl;

                        // Counter
                        counterI++;

                        if (counterI % (counterN / 10) == 0
                                || counterI == counterN) {
                            PARTONS::Partons::getInstance()->getLoggerManager()->info(
                                    "main", __func__,
                                    ElemUtils::Formatter()
                                            << "Already computed: "
                                            << 100 * counterI / counterN
                                            << "%");
                        }
                    }
                }
            }

        }

        // *******************
        // Run CFF computation
        // *******************

        if (std::string(argv[1]) == "CFF") {

            //counter
            size_t counterI = 0;
            size_t counterN = (c_nXi + 1) * (c_nTAbs + 1) * (c_nQ2 + 1);

            //loop over xi
            for (size_t iXi = 0; iXi < c_nXi; iXi++) {

                double thisXi = getValue(iXi, c_nXi, c_rangeXi, true);

                //loop over t
                for (size_t iTAbs = 0; iTAbs < c_nTAbs; iTAbs++) {

                    double thisTAbs = getValue(iTAbs, c_nTAbs, c_rangeTAbs,
                            false);

                    //loop over Q2
                    for (size_t iQ2 = 0; iQ2 < c_nQ2; iQ2++) {

                        double thisQ2 = getValue(iQ2, c_nQ2, c_rangeQ2, true);

                        // Run computation
                        PARTONS::DVCSConvolCoeffFunctionResult ccfResult =
                                pDVCSConvolCoeffFunctionService->computeForOneCCFModel(
                                        PARTONS::DVCSConvolCoeffFunctionKinematic(
                                                thisXi, -1 * thisTAbs, thisQ2,
                                                thisQ2, thisQ2), pDVCSCFFModel);

                        // Write to file
                        outputFile << thisXi << "\t" << thisTAbs << "\t"
                                << thisQ2 << "\t"
                                << ccfResult.getResult(PARTONS::GPDType::H).imag() << "\t"
                                << ccfResult.getResult(PARTONS::GPDType::H).real() << "\t"
                                << ccfResult.getResult(PARTONS::GPDType::Ht).imag() << "\t"
                                << ccfResult.getResult(PARTONS::GPDType::Ht).real() << "\t"
                                << ccfResult.getResult(PARTONS::GPDType::E).imag() << "\t"
                                << ccfResult.getResult(PARTONS::GPDType::E).real() << "\t"
                                << ccfResult.getResult(PARTONS::GPDType::Et).imag() << "\t"
                                << ccfResult.getResult(PARTONS::GPDType::Et).real()
                                << std::endl;

                        // Counter
                        counterI++;

                        if (counterI % (counterN / 10) == 0
                                || counterI == counterN) {
                            PARTONS::Partons::getInstance()->getLoggerManager()->info(
                                    "main", __func__,
                                    ElemUtils::Formatter()
                                            << "Already computed: "
                                            << 100 * counterI / counterN
                                            << "%");
                        }
                    }
                }
            }
        }

        // *******************
        // Run observable computation asymmetry
        // *******************

        if (std::string(argv[1]) == "OBS_ALU") {

            //counter
            size_t counterI = 0;
            size_t counterN = (c_nXb + 1) * (c_nTAbs + 1) * (c_nPhi + 1);

            //loop over xB
            for (size_t iXb = 0; iXb <= c_nXb; iXb++) {

                double thisXb = getValue(iXb, c_nXb, c_rangeXb, true);

                //fixed value t
                for (size_t iTAbs = 0; iTAbs <= c_nTAbs; iTAbs++) {

                    double thisTAbs = getValue(iTAbs, c_nTAbs, c_rangeTAbs,
                            false);

                    //check if valid
                    if (!isDVCSKinematicsValid(thisXb, -1 * thisTAbs, c_Q2,
                            c_E))
                        continue;

                    //loop over phi
                    for (size_t iPhi = 0; iPhi <= c_nPhi; iPhi++) {

                        double thisPhi = getValue(iPhi, c_nPhi, c_rangePhi,
                                false);

                        // Run computation
                        PARTONS::ObservableResult obsResult =
                                pObservableService->computeObservable(
                                        PARTONS::ObservableKinematic(thisXb,
                                                -1 * thisTAbs, c_Q2, c_E,
                                                thisPhi), pObservable[0]);

                        // Write to file
                        outputFile << thisXb << "\t" << thisTAbs << "\t"
                                << thisPhi << "\t" << obsResult.getValue()
                                << std::endl;

                        // Counter
                        counterI++;

                        if (counterI % (counterN / 10) == 0
                                || counterI == counterN) {
                            PARTONS::Partons::getInstance()->getLoggerManager()->info(
                                    "main", __func__,
                                    ElemUtils::Formatter()
                                            << "Already computed: "
                                            << 100 * counterI / counterN
                                            << "%");
                        }
                    }
                }
            }
        }

        // *******************
        // Run observable computation cross section
        // *******************

        if (std::string(argv[1]) == "OBS_CS") {

            //counter
            size_t counterI = 0;
            size_t counterN = (c_nXb + 1) * (c_TAbs + 1) * 3;

            //loop over xB
            for (size_t iXb = 0; iXb <= c_nXb; iXb++) {

                double thisXb = getValue(iXb, c_nXb, c_rangeXb, true);

                //check if valid
                if (!isDVCSKinematicsValid(thisXb, -1 * c_TAbs, c_Q2, c_E))
                    continue;

                //loop over phi
                for (size_t iPhi = 0; iPhi <= c_nPhi; iPhi++) {

                    double thisPhi = getValue(iPhi, c_nPhi, c_rangePhi, false);

                    PARTONS::ObservableResult obsResult[3];

                    for (size_t i = 0; i < 3; i++) {

                        // Run computation
                        obsResult[i] = pObservableService->computeObservable(
                                PARTONS::ObservableKinematic(thisXb,
                                        -1 * c_TAbs, c_Q2, c_E, thisPhi),
                                pObservable[i]);

                        // Write to file
                        outputFile << thisXb << "\t" << c_TAbs << "\t" << i
                                << "\t" << obsResult[i].getValue() << std::endl;

                        // Counter
                        counterI++;

                        if (counterI % (counterN / 10) == 0
                                || counterI == counterN) {
                            PARTONS::Partons::getInstance()->getLoggerManager()->info(
                                    "main", __func__,
                                    ElemUtils::Formatter()
                                            << "Already computed: "
                                            << 100 * counterI / counterN
                                            << "%");
                        }
                    }
                }
            }
        }

        PARTONS::Partons::getInstance()->getLoggerManager()->info("main",
                __func__,
                ElemUtils::Formatter() << "Already computed: Finished!");

        // Close file
        outputFile.flush();
        outputFile.close();

        // Remove pointer references
        // Module pointers are managed by PARTONS
        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pGPDModule, 0);
        pGPDModule = 0;

        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pDVCSCFFModel, 0);
        pDVCSCFFModel = 0;

        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pXiConverterModule, 0);
        pXiConverterModule = 0;

        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pScalesModule, 0);
        pScalesModule = 0;

        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pProcessModule, 0);
        pProcessModule = 0;

        for (size_t i = 0; i < 3; i++) {

            if (pObservable[i] != 0) {
                PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                        pObservable[i], 0);
                pObservable[i] = 0;
            }
        }

    }
// Appropriate catching of exceptions is crucial for working of PARTONS.
// PARTONS defines its own type of exception, which allows to display class name and function name
// where the exception has occurred, but also a human readable explanation.
    catch (const ElemUtils::CustomException &e) {

        // Display what happened
        pPartons->getLoggerManager()->error(e);

        // Close PARTONS application properly
        if (pPartons) {
            pPartons->close();
        }
    }
// In a case of standard exception.
    catch (const std::exception &e) {

        // Display what happened
        pPartons->getLoggerManager()->error("main", __func__, e.what());

        // Close PARTONS application properly
        if (pPartons) {
            pPartons->close();
        }
    }

// Close PARTONS application properly
    if (pPartons) {
        pPartons->close();
    }

    return 0;
}

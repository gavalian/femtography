#include "../include/examples.h"

#include <ElementaryUtils/logger/LoggerManager.h>
#include <ElementaryUtils/parameters/Parameter.h>
#include <ElementaryUtils/parameters/Parameters.h>
#include <NumA/integration/one_dimension/IntegratorType1D.h>
#include <NumA/integration/one_dimension/QuadratureIntegrator1D.h>
#include <partons/beans/convol_coeff_function/DVCS/DVCSConvolCoeffFunctionKinematic.h>
#include <partons/beans/gpd/GPDKinematic.h>
#include <partons/beans/KinematicUtils.h>
#include <partons/beans/List.h>
#include <partons/beans/observable/ObservableKinematic.h>
#include <partons/beans/PerturbativeQCDOrderType.h>
#include <partons/modules/active_flavors_thresholds/ActiveFlavorsThresholdsConstant.h>
#include <partons/modules/convol_coeff_function/DVCS/DVCSCFFStandard.h>
#include <partons/modules/evolution/gpd/GPDEvolutionVinnikov.h>
#include <partons/modules/gpd/GPDGK16Numerical.h>
#include <partons/modules/observable/DVCS/asymmetry/DVCSAllMinus.h>
#include <partons/modules/process/DVCS/DVCSProcessBMJ12.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongVinnikov.h>
#include <partons/modules/scales/ScalesQ2Multiplier.h>
#include <partons/modules/xi_converter/XiConverterXBToXi.h>
#include <partons/ModuleObjectFactory.h>
#include <partons/Partons.h>
#include <partons/services/ConvolCoeffFunctionService.h>
#include <partons/services/GPDService.h>
#include <partons/services/ObservableService.h>
#include <partons/ServiceObjectRegistry.h>

void computeSingleKinematicsForGPD() {

    // Retrieve GPD service
    PARTONS::GPDService* pGPDService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

    // Create GPD module with the BaseModuleFactory
    PARTONS::GPDModule* pGPDModel =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16Numerical::classId);

    // Create a GPDKinematic(x, xi, t, MuF2, MuR2) to compute
    PARTONS::GPDKinematic gpdKinematic(0.1, 0.2, -0.1, 2., 2.);

    // Run computation
    PARTONS::GPDResult gpdResult = pGPDService->computeGPDModel(gpdKinematic,
            pGPDModel);

    // Print results
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            gpdResult.toString());

    // Remove pointer references
    // Module pointers are managed by PARTONS
    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pGPDModel, 0);
    pGPDModel = 0;
}

void computeManyKinematicsForGPD() {

    // Retrieve GPD service
    PARTONS::GPDService* pGPDService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

    // Create GPD module with the BaseModuleFactory
    PARTONS::GPDModule* pGPDModel =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16Numerical::classId);

    // Load list of kinematics from file
    PARTONS::List<PARTONS::GPDKinematic> gpdKinematicList =
            PARTONS::KinematicUtils().getGPDKinematicFromFile(
                    "/root/workspace/partons-example/data/examples/gpd/kinematics_gpd.csv");

    // Run computation
    PARTONS::List<PARTONS::GPDResult> gpdResultList =
            pGPDService->computeManyKinematicOneModel(gpdKinematicList,
                    pGPDModel);

    // Print results
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            gpdResultList.toString());

    // Remove pointer references
    // Module pointers are managed by PARTONS
    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pGPDModel, 0);
    pGPDModel = 0;
}

void computeSingleKinematicsForDVCSComptonFormFactor() {

    // Retrieve service
    PARTONS::ConvolCoeffFunctionService* pDVCSConvolCoeffFunctionService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getConvolCoeffFunctionService();

    // Create GPD module with the BaseModuleFactory
    PARTONS::GPDModule* pGPDModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16Numerical::classId);

    // Create CFF module with the BaseModuleFactory
    PARTONS::DVCSConvolCoeffFunctionModule* pDVCSCFFModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSConvolCoeffFunctionModule(
                    PARTONS::DVCSCFFStandard::classId);

    // Create parameters to configure later DVCSCFFModel with PerturbativeQCD = LO
    ElemUtils::Parameters parameters(
            PARTONS::PerturbativeQCDOrderType::PARAMETER_NAME_PERTURBATIVE_QCD_ORDER_TYPE,
            PARTONS::PerturbativeQCDOrderType::LO);

    // Configure DVCSCFFModule with previous parameters.
    pDVCSCFFModule->configure(parameters);

    // Link modules (set physics assumptions of your computation)
    pDVCSCFFModule->setGPDModule(pGPDModule);

    // Create kinematic
    PARTONS::DVCSConvolCoeffFunctionKinematic cffKinematic =
            PARTONS::DVCSConvolCoeffFunctionKinematic(0.01, -0.1, 4., 4., 4.);

    // Run computation
    PARTONS::DVCSConvolCoeffFunctionResult cffResult =
            pDVCSConvolCoeffFunctionService->computeForOneCCFModel(cffKinematic,
                    pDVCSCFFModule);

    // Print results for DVCSCFFModule
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            cffResult.toString());

    // Remove pointer references
    // Module pointers are managed by PARTONS
    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pDVCSCFFModule, 0);
    pDVCSCFFModule = 0;

    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pGPDModule, 0);
    pGPDModule = 0;
}

void computeManyKinematicsForDVCSComptonFormFactor() {

    // Retrieve service
    PARTONS::ConvolCoeffFunctionService* pDVCSConvolCoeffFunctionService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getConvolCoeffFunctionService();

    // Create GPD module with the BaseModuleFactory
    PARTONS::GPDModule* pGPDModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16Numerical::classId);

    // Create CFF module with the BaseModuleFactory
    PARTONS::DVCSConvolCoeffFunctionModule* pDVCSCFFModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSConvolCoeffFunctionModule(
                    PARTONS::DVCSCFFStandard::classId);

    // Create parameters to configure later DVCSCFFModel with PerturbativeQCD = LO
    ElemUtils::Parameters parameters(
            PARTONS::PerturbativeQCDOrderType::PARAMETER_NAME_PERTURBATIVE_QCD_ORDER_TYPE,
            PARTONS::PerturbativeQCDOrderType::LO);

    // Configure DVCSCFFModule with previous parameters.
    pDVCSCFFModule->configure(parameters);

    // Link modules (set physics assumptions of your computation)
    pDVCSCFFModule->setGPDModule(pGPDModule);

    // Load list of kinematics from file
    PARTONS::List<PARTONS::DVCSConvolCoeffFunctionKinematic> cffKinematicList =
            PARTONS::KinematicUtils().getCCFKinematicFromFile(
                    "/root/workspace/partons-example/data/examples/cff/kinematics_dvcs_cff.csv");

    // Run computation
    PARTONS::List<PARTONS::DVCSConvolCoeffFunctionResult> cffResultList =
            pDVCSConvolCoeffFunctionService->computeForOneCCFModelAndManyKinematics(
                    cffKinematicList, pDVCSCFFModule);

    // Print results for DVCSCFFModule
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            cffResultList.toString());

    // Remove pointer references
    // Module pointers are managed by PARTONS
    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pDVCSCFFModule, 0);
    pDVCSCFFModule = 0;

    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pGPDModule, 0);
    pGPDModule = 0;
}

void computeSingleKinematicsForDVCSObservable() {

    // Retrieve Observable service
    PARTONS::ObservableService* pObservableService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getObservableService();

    // Create GPDModule
    PARTONS::GPDModule* pGPDModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16Numerical::classId);

    // Create CFF module
    PARTONS::DVCSConvolCoeffFunctionModule* pDVCSCFFModel =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSConvolCoeffFunctionModule(
                    PARTONS::DVCSCFFStandard::classId);

    // Set its PerturbativeQCDOrder
    pDVCSCFFModel->configure(
            ElemUtils::Parameter(
                    PARTONS::PerturbativeQCDOrderType::PARAMETER_NAME_PERTURBATIVE_QCD_ORDER_TYPE,
                    PARTONS::PerturbativeQCDOrderType::NLO));

    // Create XiConverterModule
    PARTONS::XiConverterModule* pXiConverterModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newXiConverterModule(
                    PARTONS::XiConverterXBToXi::classId);

    // Create ScalesModule
    PARTONS::ScalesModule* pScalesModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newScalesModule(
                    PARTONS::ScalesQ2Multiplier::classId);

    // Set its lambda parameter, so MuF2 = MuR2 = lambda * Q2
    pScalesModule->configure(
            ElemUtils::Parameter(
                    PARTONS::ScalesQ2Multiplier::PARAMETER_NAME_LAMBDA, 1.));

    // Create ProcessModule
    PARTONS::DVCSProcessModule* pProcessModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSProcessModule(
                    PARTONS::DVCSProcessBMJ12::classId);

    // Create Observable
    PARTONS::Observable* pObservable =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newObservable(
                    PARTONS::DVCSAllMinus::classId);

    // Link modules (set physics assumptions of your computation)
    pObservable->setProcessModule(pProcessModule);
    pProcessModule->setScaleModule(pScalesModule);
    pProcessModule->setXiConverterModule(pXiConverterModule);
    pProcessModule->setConvolCoeffFunctionModule(pDVCSCFFModel);
    pDVCSCFFModel->setGPDModule(pGPDModule);

    // Load list of kinematics from file
    PARTONS::ObservableKinematic observableKinematic =
            PARTONS::ObservableKinematic(0.2, -0.1, 2., 6.);

    // Create kinematic
    PARTONS::ObservableResult observableResult =
            pObservableService->computeObservable(observableKinematic,
                    pObservable);

    // Print results
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            observableResult.toString());

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

    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pObservable, 0);
    pObservable = 0;
}

void computeManyKinematicsForDVCSObservable() {

    // Retrieve Observable service
    PARTONS::ObservableService* pObservableService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getObservableService();

    // Create GPDModule
    PARTONS::GPDModule* pGPDModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16Numerical::classId);

    // Create CFF module
    PARTONS::DVCSConvolCoeffFunctionModule* pDVCSCFFModel =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSConvolCoeffFunctionModule(
                    PARTONS::DVCSCFFStandard::classId);

    // Set its PerturbativeQCDOrder
    pDVCSCFFModel->configure(
            ElemUtils::Parameter(
                    PARTONS::PerturbativeQCDOrderType::PARAMETER_NAME_PERTURBATIVE_QCD_ORDER_TYPE,
                    PARTONS::PerturbativeQCDOrderType::NLO));

    // Create XiConverterModule
    PARTONS::XiConverterModule* pXiConverterModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newXiConverterModule(
                    PARTONS::XiConverterXBToXi::classId);

    // Create ScalesModule
    PARTONS::ScalesModule* pScalesModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newScalesModule(
                    PARTONS::ScalesQ2Multiplier::classId);

    // Set its lambda parameter, so MuF2 = MuR2 = lambda * Q2
    pScalesModule->configure(
            ElemUtils::Parameter(
                    PARTONS::ScalesQ2Multiplier::PARAMETER_NAME_LAMBDA, 1.));

    // Create ProcessModule
    PARTONS::DVCSProcessModule* pProcessModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSProcessModule(
                    PARTONS::DVCSProcessBMJ12::classId);

    // Create Observable
    PARTONS::Observable* pObservable =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newObservable(
                    PARTONS::DVCSAllMinus::classId);

    // Link modules (set physics assumptions of your computation)
    pObservable->setProcessModule(pProcessModule);
    pProcessModule->setScaleModule(pScalesModule);
    pProcessModule->setXiConverterModule(pXiConverterModule);
    pProcessModule->setConvolCoeffFunctionModule(pDVCSCFFModel);
    pDVCSCFFModel->setGPDModule(pGPDModule);

    // Load list of kinematics from file
    PARTONS::List<PARTONS::ObservableKinematic> observableKinematicList =
            PARTONS::KinematicUtils().getObservableKinematicFromFile(
                    "/home/partons/git/partons-example/data/examples/observable/kinematics_dvcs_observable.csv");

    // Run computation
    PARTONS::List<PARTONS::ObservableResult> observableResultList =
            pObservableService->computeManyKinematicOneModel(
                    observableKinematicList, pObservable);

    // Print results
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            observableResultList.toString());

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

    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pObservable, 0);
    pObservable = 0;
}

void changeIntegrationRoutine() {

    // Retrieve GPD service
    PARTONS::GPDService* pGPDService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

    // Create GPD module with the BaseModuleFactory
    PARTONS::GPDModule* pGPDModel =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16Numerical::classId);

    // Change integration routine
    pGPDModel->configure(
            ElemUtils::Parameter(
                    PARTONS::MathIntegratorModule::PARAM_NAME_INTEGRATOR_TYPE,
                    NumA::IntegratorType1D::GK21_ADAPTIVE));

    // Create a GPDKinematic(x, xi, t, MuF2, MuR2) to compute
    PARTONS::GPDKinematic gpdKinematic(0.1, 0.2, -0.1, 2., 2.);

    // Run computation
    PARTONS::GPDResult gpdResult = pGPDService->computeGPDModel(gpdKinematic,
            pGPDModel);

    // Print results
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            gpdResult.toString());

    // Remove pointer references
    // Module pointers are managed by PARTONS
    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pGPDModel, 0);
    pGPDModel = 0;
}

void makeUseOfGPDEvolution() {

    // Retrieve GPD service
    PARTONS::GPDService* pGPDService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

    // Create GPD module with the BaseModuleFactory
    PARTONS::GPDModule* pGPDModel =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16Numerical::classId);

    // Create GPD evolution module with the BaseModuleFactory
    PARTONS::GPDEvolutionModule* pGPDEvolutionModel =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDEvolutionModule(
                    PARTONS::GPDEvolutionVinnikov::classId);

    // Create alphaS module with the BaseModuleFactory
    PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newRunningAlphaStrongModule(
                    PARTONS::RunningAlphaStrongVinnikov::classId);

    // Create active flavors thresholds module with the BaseModuleFactory
    PARTONS::ActiveFlavorsThresholdsModule* pActiveFlavorsThresholdsModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newActiveFlavorsThresholdsModule(
                    PARTONS::ActiveFlavorsThresholdsConstant::classId);

    // ActiveFlavorsThresholdsConstant module allows you to set the desired nf value that will be constant during performing the evolution.
    // Default value is nf = 3. You can change it in the following way, but you must be sure that both used GPD model and evolution routine can handle it.
    static_cast<PARTONS::ActiveFlavorsThresholdsConstant*>(pActiveFlavorsThresholdsModule)->setNFlavors(
            3);

    // Create parameters to configure later GPDEvolutionModule
    ElemUtils::Parameters parameters;

    // Number of steps in the factorization scale (i.e. set the number of steps in the integration over factorization scale)
    // One step is a typical value for Vinnikov code
    parameters.add(NumA::QuadratureIntegrator1D::PARAM_NAME_N, 1);

    // PerturbativeQCD = LO
    parameters.add(
            PARTONS::PerturbativeQCDOrderType::PARAMETER_NAME_PERTURBATIVE_QCD_ORDER_TYPE,
            PARTONS::PerturbativeQCDOrderType::LO);

    // Configure GPDEvolutionModule with previous parameters.
    pGPDEvolutionModel->configure(parameters);

    // Link modules (set physics assumptions of your computation)
    pGPDEvolutionModel->setRunningAlphaStrongModule(pRunningAlphaStrongModule);
    pGPDEvolutionModel->setActiveFlavorsModule(pActiveFlavorsThresholdsModule);
    pGPDModel->setEvolQcdModule(pGPDEvolutionModel);

    // Create a GPDKinematic(x, xi, t, MuF2, MuR2) to compute
    PARTONS::GPDKinematic gpdKinematic(0.1, 0.2, -0.1, 40., 40.);

    // Run computation
    PARTONS::GPDResult gpdResult = pGPDService->computeGPDModel(gpdKinematic,
            pGPDModel);

    // Print results
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            gpdResult.toString());

    // Remove pointer references
    // Module pointers are managed by PARTONS
    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pActiveFlavorsThresholdsModule, 0);
    pGPDModel = 0;

    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pRunningAlphaStrongModule, 0);
    pGPDModel = 0;

    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pGPDEvolutionModel, 0);
    pGPDModel = 0;

    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pGPDModel, 0);
    pGPDModel = 0;
}

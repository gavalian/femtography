#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <partons/Partons.h>
#include <partons/services/automation/AutomationService.h>
#include <partons/ServiceObjectRegistry.h>
#include <QtCore/qcoreapplication.h>
#include <string>
#include <vector>

#include "../include/examples.h"

/*
 * Parse XML scenarios.
 */
std::vector<std::string> parseArguments(int argc, char** argv) {
  std::vector<std::string> xmlScenarios(argc - 1);

  for (unsigned int i = 1; i < argc; i++) {
    xmlScenarios[i - 1] = argv[i];
  }

  return xmlScenarios;
}

/*
 * Main function.
 */

int main(int argc, char** argv) {

  // Init Qt4
  QCoreApplication a(argc, argv);
  PARTONS::Partons* pPartons = 0;

  try {

    // Init PARTONS application
    pPartons = PARTONS::Partons::getInstance();
    pPartons->init(argc, argv);

    // ******************************************************
    // RUN XML SCENARIO *************************************
    // ******************************************************

    // You need to provide at least one scenario via executable argument
    /*if (argc <= 1) {

      throw ElemUtils::CustomException("main", __func__,
				       "Missing argument, please provide one or more than one XML scenario file.");
    }

    // Parse arguments to retrieve xml file path list.
    std::vector<std::string> xmlScenarioFilePathList = parseArguments(argc,
								      argv);

    // Retrieve automation service parse scenario xml file and play it.
    PARTONS::AutomationService* pAutomationService =
      pPartons->getServiceObjectRegistry()->getAutomationService();

    for (unsigned int i = 0; i < xmlScenarioFilePathList.size(); i++) {
      PARTONS::Scenario* pScenario = pAutomationService->parseXMLFile(
								      xmlScenarioFilePathList[i]);
      pAutomationService->playScenario(pScenario);
      }*/
    // ******************************************************
    // RUN CPP CODE *****************************************
    // ******************************************************

    // You can put your own code here and build a stand-alone program based on PARTONS library.
    // To learn how you can use PARTONS library study provided examples of functions to be found in
    // include/examples.h (header) and src/examples.cpp (source) files.
    // To run these examples just call them here, e.g.:

    computeManyKinematicsForDVCSComptonFormFactor();
    //computeManyKinematicsForGPD();
    //computeSingleKinematicsForGPD();

    // Note, that you may need to comment out the part responsible for the running of XML scenarios.

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


#include "../include/DVCSCrossSectionUUMinusBHSubProc.h"

#include <NumA/linear_algebra/vector/Vector3D.h>
#include <partons/beans/observable/ObservableChannel.h>
#include <partons/beans/process/DVCSSubProcessType.h>
#include <partons/BaseObjectRegistry.h>
#include <partons/FundamentalPhysicalConstants.h>
#include <partons/modules/process/DVCS/DVCSProcessModule.h>

namespace PARTONS {

const unsigned int DVCSCrossSectionUUMinusBHSubProc::classId =
        BaseObjectRegistry::getInstance()->registerBaseObject(
                new DVCSCrossSectionUUMinusBHSubProc("DVCSCrossSectionUUMinusBHSubProc"));

DVCSCrossSectionUUMinusBHSubProc::DVCSCrossSectionUUMinusBHSubProc(
        const std::string &className) :
        Observable(className) {
    m_channel = ObservableChannel::DVCS;
}

DVCSCrossSectionUUMinusBHSubProc::DVCSCrossSectionUUMinusBHSubProc(
        const DVCSCrossSectionUUMinusBHSubProc& other) :
        Observable(other) {
}

DVCSCrossSectionUUMinusBHSubProc::~DVCSCrossSectionUUMinusBHSubProc() {
}

DVCSCrossSectionUUMinusBHSubProc* DVCSCrossSectionUUMinusBHSubProc::clone() const {
    return new DVCSCrossSectionUUMinusBHSubProc(*this);
}

double DVCSCrossSectionUUMinusBHSubProc::computePhiObservable(double phi) {

    double result = 0.;

    double A = static_cast<DVCSProcessModule*>(m_pProcessModule)->computeCrossSection(1., -1,
            NumA::Vector3D(0., 0., 0.), phi,
            DVCSSubProcessType::BH);

    double B = static_cast<DVCSProcessModule*>(m_pProcessModule)->computeCrossSection(-1., -1,
            NumA::Vector3D(0., 0., 0.), phi,
            DVCSSubProcessType::BH);

    result = (A + B) / 2.;

    //integrate over transversely polarized target dependence to obtain 4-fold differential cross-section
    result *= 2 * Constant::PI;

    //change to nb
    result *= Constant::CONV_GEVm2_TO_NBARN;

    return result;
}

} /* namespace PARTONS */

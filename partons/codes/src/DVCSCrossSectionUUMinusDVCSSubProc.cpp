#include "../include/DVCSCrossSectionUUMinusDVCSSubProc.h"

#include <NumA/linear_algebra/vector/Vector3D.h>
#include <partons/beans/observable/ObservableChannel.h>
#include <partons/beans/process/DVCSSubProcessType.h>
#include <partons/BaseObjectRegistry.h>
#include <partons/FundamentalPhysicalConstants.h>
#include <partons/modules/process/DVCS/DVCSProcessModule.h>

namespace PARTONS {

const unsigned int DVCSCrossSectionUUMinusDVCSSubProc::classId =
        BaseObjectRegistry::getInstance()->registerBaseObject(
                new DVCSCrossSectionUUMinusDVCSSubProc("DVCSCrossSectionUUMinusDVCSSubProc"));

DVCSCrossSectionUUMinusDVCSSubProc::DVCSCrossSectionUUMinusDVCSSubProc(
        const std::string &className) :
        Observable(className) {
    m_channel = ObservableChannel::DVCS;
}

DVCSCrossSectionUUMinusDVCSSubProc::DVCSCrossSectionUUMinusDVCSSubProc(
        const DVCSCrossSectionUUMinusDVCSSubProc& other) :
        Observable(other) {
}

DVCSCrossSectionUUMinusDVCSSubProc::~DVCSCrossSectionUUMinusDVCSSubProc() {
}

DVCSCrossSectionUUMinusDVCSSubProc* DVCSCrossSectionUUMinusDVCSSubProc::clone() const {
    return new DVCSCrossSectionUUMinusDVCSSubProc(*this);
}

double DVCSCrossSectionUUMinusDVCSSubProc::computePhiObservable(double phi) {

    double result = 0.;

    double A = static_cast<DVCSProcessModule*>(m_pProcessModule)->computeCrossSection(1., -1,
            NumA::Vector3D(0., 0., 0.), phi,
            DVCSSubProcessType::DVCS);

    double B = static_cast<DVCSProcessModule*>(m_pProcessModule)->computeCrossSection(-1., -1,
            NumA::Vector3D(0., 0., 0.), phi,
            DVCSSubProcessType::DVCS);

    result = (A + B) / 2.;

    //integrate over transversely polarized target dependence to obtain 4-fold differential cross-section
    result *= 2 * Constant::PI;

    //change to nb
    result *= Constant::CONV_GEVm2_TO_NBARN;

    return result;
}

} /* namespace PARTONS */

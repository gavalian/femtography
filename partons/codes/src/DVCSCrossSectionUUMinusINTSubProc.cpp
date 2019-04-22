#include "../include/DVCSCrossSectionUUMinusINTSubProc.h"

#include <NumA/linear_algebra/vector/Vector3D.h>
#include <partons/beans/observable/ObservableChannel.h>
#include <partons/beans/process/DVCSSubProcessType.h>
#include <partons/BaseObjectRegistry.h>
#include <partons/FundamentalPhysicalConstants.h>
#include <partons/modules/process/DVCS/DVCSProcessModule.h>

namespace PARTONS {

const unsigned int DVCSCrossSectionUUMinusINTSubProc::classId =
        BaseObjectRegistry::getInstance()->registerBaseObject(
                new DVCSCrossSectionUUMinusINTSubProc("DVCSCrossSectionUUMinusINTSubProc"));

DVCSCrossSectionUUMinusINTSubProc::DVCSCrossSectionUUMinusINTSubProc(
        const std::string &className) :
        Observable(className) {
    m_channel = ObservableChannel::DVCS;
}

DVCSCrossSectionUUMinusINTSubProc::DVCSCrossSectionUUMinusINTSubProc(
        const DVCSCrossSectionUUMinusINTSubProc& other) :
        Observable(other) {
}

DVCSCrossSectionUUMinusINTSubProc::~DVCSCrossSectionUUMinusINTSubProc() {
}

DVCSCrossSectionUUMinusINTSubProc* DVCSCrossSectionUUMinusINTSubProc::clone() const {
    return new DVCSCrossSectionUUMinusINTSubProc(*this);
}

double DVCSCrossSectionUUMinusINTSubProc::computePhiObservable(double phi) {

    double result = 0.;

    double A = static_cast<DVCSProcessModule*>(m_pProcessModule)->computeCrossSection(1., -1,
            NumA::Vector3D(0., 0., 0.), phi,
            DVCSSubProcessType::INT);

    double B = static_cast<DVCSProcessModule*>(m_pProcessModule)->computeCrossSection(-1., -1,
            NumA::Vector3D(0., 0., 0.), phi,
            DVCSSubProcessType::INT);

    result = (A + B) / 2.;

    //integrate over transversely polarized target dependence to obtain 4-fold differential cross-section
    result *= 2 * Constant::PI;

    //change to nb
    result *= Constant::CONV_GEVm2_TO_NBARN;

    return result;
}

} /* namespace PARTONS */

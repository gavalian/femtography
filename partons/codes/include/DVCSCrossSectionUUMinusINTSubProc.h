#ifndef DVCSCrossSectionUUMinusINTSubProc_H
#define DVCSCrossSectionUUMinusINTSubProc_H

#include <partons/modules/observable/Observable.h>
#include <string>

namespace PARTONS {

class DVCSCrossSectionUUMinusINTSubProc: public Observable {

public:

    /**
     * Unique ID to automatically register the class in the registry.
     */
    static const unsigned int classId;

    /**
     * Constructor.
     * @param className Name of class.
     */
    DVCSCrossSectionUUMinusINTSubProc(const std::string &className);

    /**
     * Destructor.
     */
    virtual ~DVCSCrossSectionUUMinusINTSubProc();

    virtual DVCSCrossSectionUUMinusINTSubProc* clone() const;
    virtual double computePhiObservable(double phi);

protected:

    /**
     * Copy constructor.
     * @param other Object to be copied.
     */
    DVCSCrossSectionUUMinusINTSubProc(
            const DVCSCrossSectionUUMinusINTSubProc &other);
};

} /* namespace PARTONS */

#endif /* DVCSCrossSectionUUMinusINTSubProc_H */

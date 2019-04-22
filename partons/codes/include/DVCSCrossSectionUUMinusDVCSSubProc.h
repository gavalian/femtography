#ifndef DVCSCrossSectionUUMinusDVCSSubProc_H
#define DVCSCrossSectionUUMinusDVCSSubProc_H

#include <partons/modules/observable/Observable.h>
#include <string>

namespace PARTONS {

class DVCSCrossSectionUUMinusDVCSSubProc: public Observable {

public:

    /**
     * Unique ID to automatically register the class in the registry.
     */
    static const unsigned int classId;

    /**
     * Constructor.
     * @param className Name of class.
     */
    DVCSCrossSectionUUMinusDVCSSubProc(const std::string &className);

    /**
     * Destructor.
     */
    virtual ~DVCSCrossSectionUUMinusDVCSSubProc();

    virtual DVCSCrossSectionUUMinusDVCSSubProc* clone() const;
    virtual double computePhiObservable(double phi);

protected:

    /**
     * Copy constructor.
     * @param other Object to be copied.
     */
    DVCSCrossSectionUUMinusDVCSSubProc(
            const DVCSCrossSectionUUMinusDVCSSubProc &other);
};

} /* namespace PARTONS */

#endif /* DVCSCrossSectionUUMinusDVCSSubProc_H */

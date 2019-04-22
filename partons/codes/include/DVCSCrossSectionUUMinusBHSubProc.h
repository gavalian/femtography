#ifndef DVCSCrossSectionUUMinusBHSubProc_H
#define DVCSCrossSectionUUMinusBHSubProc_H

#include <partons/modules/observable/Observable.h>
#include <string>

namespace PARTONS {

class DVCSCrossSectionUUMinusBHSubProc: public Observable {

public:

    /**
     * Unique ID to automatically register the class in the registry.
     */
    static const unsigned int classId;

    /**
     * Constructor.
     * @param className Name of class.
     */
    DVCSCrossSectionUUMinusBHSubProc(const std::string &className);

    /**
     * Destructor.
     */
    virtual ~DVCSCrossSectionUUMinusBHSubProc();

    virtual DVCSCrossSectionUUMinusBHSubProc* clone() const;
    virtual double computePhiObservable(double phi);

protected:

    /**
     * Copy constructor.
     * @param other Object to be copied.
     */
    DVCSCrossSectionUUMinusBHSubProc(
            const DVCSCrossSectionUUMinusBHSubProc &other);
};

} /* namespace PARTONS */

#endif /* DVCSCrossSectionUUMinusBHSubProc_H */

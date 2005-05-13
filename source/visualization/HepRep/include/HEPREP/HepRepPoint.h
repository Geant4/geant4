// -*- C++ -*-
// AID-GENERATED
// =========================================================================
// This class was generated by AID - Abstract Interface Definition          
// DO NOT MODIFY, but use the org.freehep.aid.Aid utility to regenerate it. 
// =========================================================================
#ifndef HEPREP_HEPREPPOINT_H
#define HEPREP_HEPREPPOINT_H 1

// Copyright 2000-2002, FreeHEP.

#include <vector>

#include "HEPREP/HepRepAttribute.h"

namespace HEPREP {

class HepRepInstance;

/**
 * HepRepPoint interface. The HepRepMath class can be used to deal with the conversions.
 *
 * @author Mark Donszelmann
 */
class HepRepPoint : virtual public HepRepAttribute {

public: 
    /// Destructor.
    virtual ~HepRepPoint() { /* nop */; }

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @return dx-coordinate
     */
    virtual double getX() = 0;

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @return dy-coordinate
     */
    virtual double getY() = 0;

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @return dz-coordinate
     */
    virtual double getZ() = 0;

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @return rho = sqrt(dx2+dy2);
     */
    virtual double getRho() = 0;

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @return phi = atan2(dy, dx);
     */
    virtual double getPhi() = 0;

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @return theta = atan2(rho, dx);
     */
    virtual double getTheta() = 0;

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @return r = sqrt(dx2+dy2+dz2);
     */
    virtual double getR() = 0;

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @return eta = -0.5*clog((1.-ct)/(1.+ct)), where ct = .cos(getTheta(dx, dy, dz));
     */
    virtual double getEta() = 0;

    /**
     * Returns coordinate with respect to vertex at (xVertex, yVertex, zVertex).
     *
     * @return dx-coordinate
     */
    virtual double getX(double xVertex, double yVertex, double zVertex) = 0;

    /**
     * Returns coordinate with respect to vertex at (xVertex, yVertex, zVertex).
     *
     * @return dy-coordinate
     */
    virtual double getY(double xVertex, double yVertex, double zVertex) = 0;

    /**
     * Returns coordinate with respect to vertex at (xVertex, yVertex, zVertex).
     *
     * @return dz-coordinate
     */
    virtual double getZ(double xVertex, double yVertex, double zVertex) = 0;

    /**
     * Returns coordinate with respect to vertex at (xVertex, yVertex, zVertex).
     *
     * @return rho = sqrt(dx2+dy2);
     */
    virtual double getRho(double xVertex, double yVertex, double zVertex) = 0;

    /**
     * Returns coordinate with respect to vertex at (xVertex, yVertex, zVertex).
     *
     * @return phi = atan2(dy, dx);
     */
    virtual double getPhi(double xVertex, double yVertex, double zVertex) = 0;

    /**
     * Returns coordinate with respect to vertex at (xVertex, yVertex, zVertex).
     *
     * @return theta = atan2(rho, dx);
     */
    virtual double getTheta(double xVertex, double yVertex, double zVertex) = 0;

    /**
     * Returns coordinate with respect to vertex at (xVertex, yVertex, zVertex).
     *
     * @return r = sqrt(dx2+dy2+dz2);
     */
    virtual double getR(double xVertex, double yVertex, double zVertex) = 0;

    /**
     * Returns coordinate with respect to vertex at (xVertex, yVertex, zVertex).
     *
     * @return eta = -0.5*clog((1.-ct)/(1.+ct)), where ct = .cos(getTheta(dx, dy, dz));
     */
    virtual double getEta(double xVertex, double yVertex, double zVertex) = 0;

    /**
     * Returns coordinate with respect to vertex at (0, 0, 0).
     *
     * @param xyz list of three coordinates which are filled and returned.
     *            If null, a new list of three coordinates is allocated.
     * @return list of 3 coordinates.
     */
    virtual std::vector<double>  * getXYZ(std::vector<double>  * xyz) = 0;

    /**
     * Returns associated instance (parent).
     *
     * @return HepRepInstance.
     */
    virtual HepRepInstance * getInstance() = 0;

    /**
     * Returns a deep copy of this point.
     *
     * @param parent to add the copy to.
     * @return copy of this point.
     */
    virtual HepRepPoint * copy(HepRepInstance * parent) = 0;
}; // class
} // namespace HEPREP
#endif /* ifndef HEPREP_HEPREPPOINT_H */

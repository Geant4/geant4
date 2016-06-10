//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UEnclosingCylinder
//
// Class description:
//
//   Definition of a utility class for quickly deciding if a point
//   is clearly outside a polyhedra or polycone or deciding if
//   a trajectory is clearly going to miss those shapes.
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UEnclosingCylinder_hh
#define UEnclosingCylinder_hh

#include "UTypes.hh"
#include "UTubs.hh"

class UReduciblePolygon;

class UEnclosingCylinder
{
  public: // with description

    UEnclosingCylinder(/*const UReduciblePolygon *rz*/  double r, double lo, double hi,
                                                        bool phiIsOpen,
                                                        double startPhi, double totalPhi);
    ~UEnclosingCylinder();

    bool MustBeOutside(const UVector3& p) const;
    // Decide very rapidly if the point is outside the cylinder.
    // If one is not certain, return false.

    bool ShouldMiss(const UVector3& p, const UVector3& v) const;
    // Decide very rapidly if the trajectory is going to miss the cylinder.
    // If one is not sure, return false.

    double DistanceTo(const UVector3& p, const UVector3& v) const;

    double SafetyFromOutside(const UVector3& p) const;

  public: // without description

    void Extent(UVector3& aMin, UVector3& aMax) const;

    double radius;    // radius of our cylinder

  protected:

    double zLo, zHi;  // z extent

    bool    phiIsOpen; // true if there is a phi segment
    double  startPhi, // for isPhiOpen==true, starting of phi segment
            totalPhi; // for isPhiOpen==true, size of phi segment

    double rx1, ry1,
           dx1, dy1;
    double rx2, ry2,
           dx2, dy2;

    bool   concave; // true, if x/y Cross section is concave

    UTubs* tube;
};

#endif

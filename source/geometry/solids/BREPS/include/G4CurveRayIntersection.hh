//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// ----------------------------------------------------------------------
// Class G4CurveRayIntersection
//
// Class Description:
//
// Class capable of storing both the parametric and the non-parametric
// representation of a intersection point on a curve. It's subclassed
// from G4CurvePoint.

// Author: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef included_G4CurveRayIntersection
#define included_G4CurveRayIntersection

#include "G4CurvePoint.hh"
#include "G4Ray.hh"

class G4CurveRayIntersection : public G4CurvePoint
{

public:  // with description

  G4CurveRayIntersection();
    // Default constructor. Sets only distance to infinity.
    // Must be followed by Init!

  G4CurveRayIntersection(G4Curve& c0, const G4Ray& r0);
    // Constructor taking a curve and a ray.
 
  ~G4CurveRayIntersection();
    // Empty destructor.

  G4CurveRayIntersection(const G4CurveRayIntersection& cr);
  G4CurveRayIntersection& operator=(const G4CurveRayIntersection& cr);
    // Copy constructor and assignment operator.

  inline void Init(G4Curve& c0, const G4Ray& r0);
    // Initialises a G4CurveRayIntersection. Called by the constructor above.

  inline const G4Ray& GetRay() const;
    // Returns the ray of intersection.

  inline void Reset();
  inline void ResetPPoint(G4double u0);
  inline void Reset(const G4Point3D& p0);
  inline void Reset(G4double u0, const G4Point3D& p0);
  inline void ResetDistance(G4double d0);
  inline void Reset(G4double u0, G4double d0);
  inline void Reset(const G4Point3D& p0, G4double d0);
  inline void Reset(G4double u0, const G4Point3D& p0, G4double d0);
    // Resets point's attributes.

  inline G4double GetPPoint();
  inline const G4Point3D& GetPoint();
    // Returns point as parameter or as point in space.

  inline G4double GetDistance();
    // Returns intersection's distance.

public:  // without description

  // additional functionalities

  inline void Update(G4CurveRayIntersection& is);
  inline void UpdateWithPointOnCurve(G4CurveRayIntersection& is);

protected:

  // data

  const G4Ray* r;
  G4double d;
  G4double kCarTolerance;
  static const G4int dFlag;

};

#include "G4CurveRayIntersection.icc"

#endif

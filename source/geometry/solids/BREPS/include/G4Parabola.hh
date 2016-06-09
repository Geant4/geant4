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
// Class G4Parabola
//
// Class description:
// 
// Definition of a generic parabola.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __PARABOLICCURVE_H
#define __PARABOLICCURVE_H 

#include "G4Conic.hh"

class G4Parabola : public G4Conic
{

public:  // with description

  G4Parabola();
  virtual ~G4Parabola();
    // Constructor & destructor.

  G4Parabola(const G4Parabola& right);
  G4Parabola& operator=(const G4Parabola& right);
    // Copy constructor and assignment operator.

  G4Curve* Project(const G4Transform3D& tr = G4Transform3D::Identity);
    // Transforms and projects the curve.

  G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Returns tangent to curve at a given point, if existing.
    // The tangent is computed from the 3D point representation.

  inline G4double  GetPMax() const;
  inline G4Point3D GetPoint(G4double param) const;
  inline G4double  GetPPoint(const G4Point3D& p) const;
    // Accessors methods.

  inline G4double GetFocalDist() const;
  inline void Init(const G4Axis2Placement3D& position0, G4double focalDist0);
    // Get/Set for the geometric data.

public:  // without description

  inline G4int IntersectRay2D(const G4Ray& ray);

  // void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  // G4Parabola(G4Point3d, G4Point3d, G4double );
  // G4Point3d EvaluateByParameterValue(G4double u);
  // G4Point3d GetBoundMax();
  // G4Point3d GetBoundMin();   

protected:

  void InitBounded();

private:

  // geometric data

  G4double focalDist;

  // for the intersection
  G4Point3D F;
  G4Point3D L0;

};

#include "G4Parabola.icc"

#endif

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
// Class G4Ellipse
//
// Class description:
// 
// Definition of a generic ellipse curve.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __ELLIPTICCURVE_H
#define __ELLIPTICCURVE_H 

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4CircularCurve.hh"

class G4Ellipse : public G4Conic
{

public:  // with description

  G4Ellipse();
  virtual ~G4Ellipse();
    // Constructor & destructor.

  G4Ellipse(const G4Ellipse& right);
  G4Ellipse& operator=(const G4Ellipse& right);
    // Copy constructor and assignment operator.

  G4Curve* Project(const G4Transform3D& tr= G4Transform3D::Identity);
    // Transforms and projects the curve.

  G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Returns tangent to curve at a given point, if existing.
    // The tangent is computed from the 3D point representation.

  inline G4double  GetPMax() const;
  inline G4Point3D GetPoint(G4double param) const;
  inline G4double  GetPPoint(const G4Point3D& p) const;
    // Accessors methods.

  inline G4double GetSemiAxis1() const;
  inline G4double GetSemiAxis2() const;
  inline void Init(const G4Axis2Placement3D& position0,
	           G4double semiAxis10, G4double semiAxis20);
    // Get/Set methods for the geometric data.


public:

  inline G4int IntersectRay2D(const G4Ray& ray);
  // void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);

protected:

  void InitBounded();

private:

  // geometric data

  G4double semiAxis1;
  G4double semiAxis2;
  G4double ratioAxis2Axis1;
  
  G4Transform3D toUnitCircle;

  G4double forTangent; // -R_1^2/R_2^2
};

#include "G4Ellipse.icc"

#endif

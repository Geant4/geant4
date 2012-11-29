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
// Class G4CircularCurve
//
// Class description:
// 
// Definition of a generic circle.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __CIRCULARCURVE_H
#define __CIRCULARCURVE_H 

#include "G4Conic.hh"

class G4CircularCurve : public G4Conic
{

public: // with description

  G4CircularCurve();
  virtual ~G4CircularCurve();
    // Constructor & destructor.

  G4CircularCurve(const G4CircularCurve& right);
  G4CircularCurve& operator=(const G4CircularCurve& right);
    // Copy-constructor and assignment operator.

  virtual G4Curve* Project(const G4Transform3D& tr= G4Transform3D::Identity);
    // Project along trasformation tr.

  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Tangent computed from the 3D point representation as for all conics.
 
  virtual G4double GetPMax() const;
  virtual G4Point3D GetPoint(G4double param) const;
  virtual G4double GetPPoint(const G4Point3D& p) const;
  inline G4double GetRadius() const;
    // Accessors for the geometric data.

  void Init(const G4Axis2Placement3D& position0, G4double radius0);
    // Initialises geometric data.

public: // without description

  // virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  virtual G4int IntersectRay2D(const G4Ray& ray);
  
protected:  // with description

  virtual void InitBounded();
    // Compute bounding box.

private:

  G4double radius;
    // Geometric data.

};

#include "G4CircularCurve.icc"

#endif

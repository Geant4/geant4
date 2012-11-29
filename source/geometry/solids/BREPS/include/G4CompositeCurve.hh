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
// Class G4CompositeCurve
//
// Class description:
// 
// Definition of a generic composition of curves.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4COMPOSITCURVE_H
#define __G4COMPOSITCURVE_H

#include "G4Curve.hh"
#include "G4CurveVector.hh"
#include "G4CurveRayIntersection.hh"
#include "G4Point3DVector.hh"


class G4CompositeCurve : public G4Curve
{

public:  // with description

  G4CompositeCurve();
  virtual ~G4CompositeCurve();
    // Default constructor & destructor

  G4CompositeCurve(const G4Point3DVector& vertices);
    // Constructor creating closed polygons, given the vertices.
    // No call to Init and SetBounds is needed after calling this
    // constructor.

  virtual G4String GetEntityType() const;
    // Returns entity type identifier.

  virtual G4Curve* Project(const G4Transform3D& tr = G4Transform3D::Identity);
    // Project along trasformation tr.
	
  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Tangent computed from the 3D point representation.

  virtual G4double  GetPMax() const;
  virtual G4Point3D GetPoint(G4double param) const;
  virtual G4double  GetPPoint(const G4Point3D& p) const;
  const G4CurveVector& GetSegments() const;
    // Acccessors for the geometric data

  void Init(const G4CurveVector& segments0);
    // Initialises geometric data.
    // The object is not responsible for deleting the curves;
    // only a shallow copy of the CurveVector is made.

public:  // without description

  // virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  virtual G4int IntersectRay2D(const G4Ray& ray);

protected:  // with description
  
  virtual void InitBounded();
    // Compute bounding box.

private:

  G4CompositeCurve(const G4CompositeCurve&);
  G4CompositeCurve& operator=(const G4CompositeCurve&);
     // Private copy-constructor and assignment operator.

private:
  
  G4CurveVector segments;
  G4CurveRayIntersection lastIntersection;
    // geometric data

};

#include "G4CompositeCurve.icc"

#endif

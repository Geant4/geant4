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
// $Id: G4BoundingBox3D.hh,v 1.8 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BoundingBox3D
//
// Class description:
// 
// Definition of a generic solid's bounding box in the 3D space.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BoundingBox3D_h
#define __G4BoundingBox3D_h 1

#include "G4Ray.hh"
#include "G4Point3D.hh"
#include "G4Vector3D.hh"

class G4BoundingBox3D
{

public:  // with description

  G4BoundingBox3D();
  G4BoundingBox3D(const G4Point3D&);
  G4BoundingBox3D(const G4Point3D&, const G4Point3D&);    
  ~G4BoundingBox3D();
    // Constructors & destructor.

  G4BoundingBox3D(const G4BoundingBox3D& right);
  G4BoundingBox3D& operator=(const G4BoundingBox3D& right);
    // Copy constructor and assignment operator.

  void Init(const G4Point3D&);
  void Init(const G4Point3D&, const G4Point3D&);
  void Extend(const G4Point3D&);
    // To create/extend the bounding box

  G4Point3D GetBoxMin() const;
  G4Point3D GetBoxMax() const;
  G4double GetDistance() const;
  void SetDistance(G4double distance0);
    // Accessors.

  G4int Inside(const G4Point3D&) const;
    // Returns 1 if the point is inside and on the bbox.
    // Returns 0 if the point is outside the bbox.

public:

  G4int GetTestResult() const;
  G4int Test(const G4Ray&);

  static const G4BoundingBox3D space;

private:

  G4int BoxIntersect(const G4Point3D&, 
		     const G4Point3D&, 
		     const G4Vector3D&) const;

  G4double DistanceToIn(const G4Point3D&,
			const G4Vector3D&) const;			  
    
private:

  G4Point3D box_min;
  G4Point3D box_max;
  G4double distance;

  G4int test_result;

  G4Point3D MiddlePoint;
  G4Vector3D GeantBox;
  G4double kCarTolerance;
};

#include "G4BoundingBox3D.icc"

#endif

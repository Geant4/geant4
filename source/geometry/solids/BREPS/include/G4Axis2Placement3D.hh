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
// Class G4Axis2Placement3D
//
// Class description:
//
// Class defining an axis placed in a 3D space.
// Attributes are: a point location, an axis vector with its
// reference direction and coordinates in space.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4Placement3D_h
#define __G4Placement3D_h 1

#include "G4Point3D.hh"
#include "G4Vector3D.hh"
#include "G4Transform3D.hh"
#include "G4PointRat.hh"
#include "G4Ray.hh"

class G4Axis2Placement3D 
{

public: // with description

  G4Axis2Placement3D();
  G4Axis2Placement3D( const G4Vector3D& refDirection0 ,
		      const G4Vector3D& axis0         , 
		      const G4Point3D&  location0      );
  ~G4Axis2Placement3D();
    // Constructors & destructor
  
  G4Axis2Placement3D(const G4Axis2Placement3D&);
  G4Axis2Placement3D& operator=(const G4Axis2Placement3D&);
    // Copy constructor and assignment operator.

  inline G4bool operator==(const G4Axis2Placement3D& other) const;
    // Equality operator.

  inline void Init( const G4Vector3D& refDirection0 , 
	            const G4Vector3D& axis0         ,
	            const G4Point3D&  location0      );

  inline G4Point3D  GetLocation() const;
  inline G4Vector3D GetAxis() const;
  inline G4Vector3D GetRefDirection() const;
    // Get/Set for geometric data

  inline G4Vector3D GetPX() const;
  inline G4Vector3D GetPY() const;
  inline G4Vector3D GetPZ() const;
    // Placement coordinate axes

  inline const G4Transform3D& GetToPlacementCoordinates() const;
  inline const G4Transform3D& GetFromPlacementCoordinates() const; 
    // Transformation from/to the placement coordinate system

public:

  //inline void Project (G4ThreeVec& Coord, const G4ThreeVec& Pt2, 
  //                     const G4Plane& Pl1, const G4Plane& Pl2)
  // {
  //   Coord.X(Pt2.X()*Pl1.a + Pt2.Y()*Pl1.b + Pt2.Z()*Pl1.c - Pl1.d);
  //   Coord.Y(Pt2.X()*Pl2.a + Pt2.Y()*Pl2.b + Pt2.Z()*Pl2.c - Pl2.d);
  //   Coord.Z(0);
  // }

private:

  // geometric data
  G4Point3D location;
  G4Vector3D axis;
  G4Vector3D refDirection;  

  // placement coordinate axes
  G4Vector3D pX, pY, pZ;

  G4Transform3D toPlacementCoordinates;
  G4Transform3D fromPlacementCoordinates;

};

#include "G4Axis2Placement3D.icc"

#endif

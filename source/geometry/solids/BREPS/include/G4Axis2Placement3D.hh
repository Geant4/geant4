// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Axis2Placement3D.hh,v 1.3 2000-11-08 14:21:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

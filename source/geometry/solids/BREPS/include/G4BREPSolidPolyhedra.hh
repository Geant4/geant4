//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4BREPSolidPolyhedra.hh,v 1.10 2002-11-06 23:28:50 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidPolyhedra
//
//  Class description:
// 
//  The polygonal solid G4BREPSolidPolyhedra is a shape defined by an inner 
//  and outer polygonal surface and two planes perpendicular to the Z axis. 
//  Each polygonal surface is created by linking a series of polygons created 
//  at different planes perpendicular to the Z-axis. All these polygons all 
//  have the same number of sides (sides) and are defined at the same Z planes 
//  for both inner and outer polygonal surfaces. 
//
//  G4BREPSolidPolyhedra( const G4String& name,
//                              G4double  phi1,
//                              G4double  dphi,
//                              G4int     sides,
//                              G4int     num_z_planes,      
//                              G4double  z_start,
//                              G4double  z_values[],
//                              G4double  RMIN[],
//                              G4double  RMAX[] )

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, S.Giani, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPPOLYHEDRA
#define __G4BREPPOLYHEDRA

#include "G4BREPSolid.hh"

class G4BREPSolidPolyhedra : public G4BREPSolid
{
public:  // with description

  G4BREPSolidPolyhedra( const G4String& name,
			G4double phi1,
			G4double dphi,
			G4int    sides,
			G4int    num_z_planes,      
			G4double z_start,
			G4double z_values[],
			G4double RMIN[],
			G4double RMAX[] );
    // Constructor defining the polyhedra through its
    // polygonal surfaces and planes.

  ~G4BREPSolidPolyhedra();
    // Destructor.

  void Initialize();
    // Calculates the bounding box.

  EInside Inside(register const G4ThreeVector& Pt) const;
    // Determines if the point Pt is inside, outside or on the surface
    // of the solid.

  G4ThreeVector SurfaceNormal(const G4ThreeVector&) const;
    // Calculates the normal of the surface at a point on the surface
    // The sense of the normal depends on the sense of the surface.

  G4double DistanceToIn(const G4ThreeVector&) const;
    // Calculates the shortest distance ("safety") from a point
    // outside the solid to any boundary of this solid.
    // Return 0 if the point is already inside.

  G4double DistanceToIn(register const G4ThreeVector& Pt, 
			register const G4ThreeVector& V) const;
    // Calculates the distance from a point Pt outside the solid
    // to the solid's boundary along a specified direction vector V.
    // Note: Intersections with boundaries less than the tolerance must
    //       be ignored if the direction is away from the boundary.

  G4double DistanceToOut(register const G4ThreeVector& Pt, 
			 register const G4ThreeVector& V, 
			 const G4bool calcNorm=false, 
			 G4bool *validNorm=0, G4ThreeVector *n=0) const;
    // Calculates the distance from a point inside the solid to the solid`s
    // boundary along a specified direction vector.
    // Return 0 if the point is already outside.
    // Note: If the shortest distance to a boundary is less than the
    //       tolerance, it is ignored. This allows for a point within a
    //       tolerant boundary to leave immediately.

  G4double DistanceToOut(const G4ThreeVector&) const;
    // Calculates the shortest distance ("safety") from a point inside the
    // solid to any boundary of this solid.
    // Return 0 if the point is already outside.	

  virtual G4std::ostream& StreamInfo(G4std::ostream& os) const;
    // Streams solid contents to output stream.

public:

  G4Polyhedron* CreatePolyhedron () const;
    // Creates a G4Polyhedron

  void Reset() const;
    // Resets all distance attributes.

private:

  G4BREPSolidPolyhedra(const G4BREPSolidPolyhedra&);
  G4BREPSolidPolyhedra& operator=(const G4BREPSolidPolyhedra&);
    // Private copy constructor and assignment operator.
  
  typedef enum _SF {
    EInverse = 0,
    ENormal = 1
  } ESurfaceSense;
    
  G4Surface* CreateTrapezoidalSurface( G4double r1, G4double r2,
                                       const G4Point3D& origin, G4double zDistance,
                                       G4Vector3D& xAxis, G4double partAngle,
                                       ESurfaceSense sense );
  
  G4Surface* CreateTriangularSurface( G4double r1, G4double r2,
                                      const G4Point3D& origin, G4double zDistance,
                                      G4Vector3D& xAxis, G4double partAngle,
                                      ESurfaceSense sense );
  
  G4Surface* ComputePlanarSurface( G4double r1, G4double r2,
                                   const G4Point3D& origin, G4Vector3D& xAxis,
                                   G4int sides, G4double partAngle,
                                   ESurfaceSense sense );

  //   The following is only utilised in storing the shape parameters for
  //  use in visualising this shape.  J.A. Feb  24, 1997
  //
  // R. Chytracek, Nov 2002, Update to new IO dumping mechanism
  
  struct {
    G4double  start_angle;
    G4double  opening_angle;
    G4int     sides;
    G4int     num_z_planes;
    G4double  z_start;               
    G4double* z_values;
    G4double* RMIN;
    G4double* RMAX;
  } constructorParams;
      
};

#endif

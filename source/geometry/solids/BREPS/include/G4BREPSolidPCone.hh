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
// Class G4BREPSolidPCone
//
// Class description:
//
//  The polyconical solid G4BREPSolidPCone is a shape defined by a set of 
//  inner and outer conical or cylindrical surface sections and two planes 
//  perpendicular to the Z axis. Each conical surface is defined by its 
//  radius at two different planes perpendicular to the Z-axis. Inner and 
//  outer conical surfaces are defined using common Z planes:
//
//  G4BREPSolidPCone( const G4String& name,
//                          G4double start_angle,
//                          G4double opening_angle,		   
//                          G4int    num_z_planes,
//                          G4double z_start,		   
//                          G4double z_values[],
//                          G4double RMIN[],
//                          G4double RMAX[] )

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, S.Giani, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidPCone
#define __G4BREPSolidPCone

#include "G4BREPSolid.hh"

class G4BREPSolidPCone : public G4BREPSolid
{

 public: // with description

  G4BREPSolidPCone( const G4String& name,
		    G4double start_angle,
		    G4double opening_angle,		   
		    G4int    num_z_planes, // sections
		    G4double z_start,		   
		    G4double z_values[],
		    G4double RMIN[],
		    G4double RMAX[] );
    // Constructor defining the polycone through its
    // conical/cylindrical surfaces.

  ~G4BREPSolidPCone();
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

public:

  void Reset() const;
    // Resets all distance attributes.

  G4Polyhedron* CreatePolyhedron () const;
    // Creates a G4Polyhedron

  G4VSolid* Clone() const;
    // Returns a pointer of a dynamically allocated copy of the solid.

  std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

public:  // without description

  G4BREPSolidPCone(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4BREPSolidPCone(const G4BREPSolidPCone& rhs);
  G4BREPSolidPCone& operator=(const G4BREPSolidPCone& rhs);
    // Copy constructor and assignment operator.

private:

  void InitializePCone();

  G4Surface* ComputePlanarSurface( G4double r1, G4double r2,
                                   G4ThreeVector& origin,
                                   G4ThreeVector& planeAxis,
                                   G4ThreeVector& planeDirection,
                                   G4int surfSense );
    // For a given radius values compute a planar surface

  // The following is only utilised in storing the shape parameters for
  // use in visualising this shape.  J.A. Feb  24, 1997
  // R. Chytracek, Nov 2002, Update to new IO dumping mechanism
  
  struct G4BREPPConeParams
  {
    G4double  start_angle;
    G4double  opening_angle;
    G4int     num_z_planes;
    G4double  z_start;               
    G4double* z_values;
    G4double* RMIN;
    G4double* RMAX;
  } constructorParams;
};

#endif

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidPolyhedra.hh,v 1.1 1999-01-07 16:07:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4BREPPOLYHEDRA
#define __G4BREPPOLYHEDRA
#include "G4BREPSolid.hh"

class G4BREPSolidPolyhedra: public G4BREPSolid
{
public:
  // Constructor for Geant3 PGon shape
  G4BREPSolidPolyhedra(
			 G4String name,
			 const G4double phi1,
			 const G4double dphi,
			 const int sides,
			 const int num_z_planes,      
			 const G4double z_start,
			 const G4double z_values[],
			 const G4double RMIN[],
			 const G4double RMAX[]     
		       );
  void Initialize();
  inline void Reset() const
    {
      Active(1);
      ((G4BREPSolidPolyhedra*)this)->intersectionDistance=kInfinity;
      StartInside(0);
      for(register int a=0;a<nb_of_surfaces;a++)
	SurfaceVec[a]->Reset();
      ShortestDistance = kInfinity;
    }

  EInside Inside(register const G4ThreeVector&) const;
  G4ThreeVector SurfaceNormal(const G4ThreeVector&) const;

  G4double DistanceToIn(const G4ThreeVector&) const;
  G4double DistanceToIn(register const G4ThreeVector&, 
			register const G4ThreeVector&) const;

  G4double DistanceToOut(register const G4ThreeVector&, 
			 register const G4ThreeVector&, 
			 const G4bool calcNorm=false, 
			 G4bool *validNorm=0, G4ThreeVector *n=0) const;
  G4double DistanceToOut(const G4ThreeVector&) const;

  ~G4BREPSolidPolyhedra();
  G4Polyhedron* CreatePolyhedron () const;

private:

  //   The following is only utilised in storing the shape parameters for
  //  use in visualising this shape.  J.A. Feb  24, 1997
  //
  struct PGonParameters {
     G4double  Start_angle;
     G4double  Opening_angle;		   
     int       Sides; 
     int       Num_z_planes; 
     // G4double z_start;		   
     G4double  *Z_values;
     G4double  *Rmin;
     G4double  *Rmax;
  }  original_parameters;
};

#endif






// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidPCone.hh,v 1.3 1999-12-15 14:49:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// The polyconical solid G4BREPSolidPCone is a shape defined by a set of 
// inner and outer conical or cylindrical surface sections and two planes 
// perpendicular to the Z axis. Each conical surface is defined by its 
// radius at two different planes perpendicular to the Z-axis. Inner and 
// outer conical surfaces are defined using common Z planes. 
//

#ifndef __G4BREPSolidPCone
#define __G4BREPSolidPCone
#include "G4BREPSolid.hh"

class G4BREPSolidPCone : public G4BREPSolid
{
 public:
  G4BREPSolidPCone( G4String name,
		    const G4double  start_angle,
		    const G4double  opening_angle,		   
		    const int       num_z_planes, // sections,
		    const G4double  z_start,		   
		    const G4double  z_values[],
		    const G4double  RMIN[],
		    const G4double  RMAX[]
		  );

  inline void Reset() const
    {
      Active(1);
      ((G4BREPSolidPCone*)this)->intersectionDistance=kInfinity;
      StartInside(0);
      for(register int a=0;a<nb_of_surfaces;a++)
	SurfaceVec[a]->Reset();
      ShortestDistance = kInfinity;
    }

  void Initialize();
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

  ~G4BREPSolidPCone();
  G4Polyhedron* CreatePolyhedron () const;

private:

  //   The following is only utilised in storing the shape parameters for
  //  use in visualising this shape.  J.A. Feb  24, 1997
  //
  struct PConeParameters {
     G4double Start_angle;
     G4double Opening_angle;		   
     int      Num_z_planes; 
     // G4double z_start;		   
     G4double *Z_values;
     G4double *Rmin;
     G4double *Rmax;
  }  original_parameters;
};

#endif

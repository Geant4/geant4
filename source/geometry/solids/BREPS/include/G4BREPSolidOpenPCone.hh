// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef G4BREPSolidOpenPCone_hh
#define G4BREPSolidOpenPCone_hh

#include "G4IntersectionSolid.hh"

class G4BREPSolidOpenPCone : public G4IntersectionSolid {

public:
  
  G4BREPSolidOpenPCone ( G4String name,
                    const G4double  start_angle,
                    const G4double  opening_angle,                 
                    const int       num_z_planes, // sections,
                    const G4double  z_start,               
                    const G4double  z_values[],
                    const G4double  RMIN[],
                    const G4double  RMAX[]
                  );

  void DescribeYourselfTo (G4VGraphicsScene& scene) const;

private:
  
};


#endif

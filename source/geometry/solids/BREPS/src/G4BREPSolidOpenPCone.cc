// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "G4BREPSolidOpenPCone.hh"
#include "G4BREPSolidPCone.hh"
#include "G4Tubs.hh"
#include "G4VGraphicsScene.hh"

G4BREPSolidOpenPCone::G4BREPSolidOpenPCone  (G4String name,
                                   const G4double start_angle,
                                   const G4double opening_angle,
                                   const int      num_z_planes, // sections,
                                   const G4double z_start,                 
                                   const G4double z_values[],
                                   const G4double RMIN[],
                                   const G4double RMAX[]
                                   ) : 
                  G4IntersectionSolid ( name, 
					new G4BREPSolidPCone ( name, start_angle, opening_angle, 
							       num_z_planes,
							       z_start, z_values, RMIN, RMAX ) ,
					new G4Tubs( "IntersectionTubs", 0, 1*cm, 1*cm, 0*deg, 
						    360*deg ) )
{


// compute max radius

  G4double MaxRMAX = 0;
  for ( int i = 0; i < num_z_planes; i++ ) 
    if ( RMAX[i] > MaxRMAX ) MaxRMAX = RMAX[i];
  		
  G4double length = z_values[num_z_planes-1] - z_values[0];
  
  G4Tubs *ptrTubs = (G4Tubs*) fPtrSolidB;  // To be modified: TODO

  ptrTubs->SetInnerRadius ( 0 );
  ptrTubs->SetOuterRadius ( MaxRMAX );
  ptrTubs->SetZHalfLength ( length );
  ptrTubs->SetStartPhiAngle ( start_angle );
  ptrTubs->SetDeltaPhiAngle ( opening_angle );

}

#include "G4Polyhedron.hh"   


void G4BREPSolidOpenPCone::DescribeYourselfTo (G4VGraphicsScene& scene) const {
  scene.AddThis ( *fPtrSolidA );
}

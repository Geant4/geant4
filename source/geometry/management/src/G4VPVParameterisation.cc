// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPVParameterisation.cc,v 1.1 1999-01-07 16:07:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Default implementations for Parameterisations that do not
//    parameterise solid and/or material
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
// #include "G4Material.hh"

G4VSolid* G4VPVParameterisation::ComputeSolid(const G4int,
                                       G4VPhysicalVolume *pPhysicalVol) 
{
    return pPhysicalVol->GetLogicalVolume()->GetSolid();
}
				       
G4Material* G4VPVParameterisation::ComputeMaterial(const G4int,
                                       G4VPhysicalVolume *pPhysicalVol) 
{
    return pPhysicalVol->GetLogicalVolume()->GetMaterial();
}				       

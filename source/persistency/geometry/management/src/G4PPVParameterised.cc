// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPVParameterised.cc,v 1.2 1999/12/15 14:51:24 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
//
//

#include "G4PPVParameterised.hh"

#include "HepODBMS/clustering/HepClustering.h"

#include "globals.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4PVParameterised.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PLogicalVolume.hh"

G4PPVParameterised::G4PPVParameterised( G4VPhysicalVolume* PhysVol,
		                    HepRef(G4PLogicalVolume) persLogVol)
 :  G4PPVReplica(PhysVol, persLogVol)
{
  // G4VPVParameterisation class has no data member
  //  fparam = NULL;
}

G4PPVParameterised::~G4PPVParameterised()
{
}

G4VPhysicalVolume* G4PPVParameterised::MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother)
{
  const G4String pName = (const char *) GetName();
  EAxis pAxis = faxis;
  G4int nReplicas = fnReplicas;

  G4VPhysicalVolume* aPhysVol = new G4PVParameterised(
       pName, aLogical, aMother, pAxis, nReplicas, NULL);

  return aPhysVol;
}

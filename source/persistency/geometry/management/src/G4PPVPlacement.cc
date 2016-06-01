// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPVPlacement.cc,v 2.5 1998/11/10 18:31:06 morita Exp $
// GEANT4 tag $Name: geant4-00 $
//
//
//

#include "G4PPVPlacement.hh"

#include "HepODBMS/clustering/HepClustering.h"

#include "globals.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PLogicalVolume.hh"

G4PPVPlacement::G4PPVPlacement(
                 G4VPhysicalVolume* PhysVol,
		         HepRef(G4PLogicalVolume) persLogVol)
 :  G4PVPhysicalVolume(PhysVol, persLogVol)
{
  fcopyNo = PhysVol->GetCopyNo();
  fmany   = PhysVol->IsMany();  
}

G4PPVPlacement::~G4PPVPlacement()
{
}

G4VPhysicalVolume* G4PPVPlacement::MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother )
{
  G4RotationMatrix* pRot = GetRotation();
  G4ThreeVector tlate = GetTranslation();
  G4String pName;

  G4VPhysicalVolume* aPhysVol = new G4PVPlacement(
                pRot, tlate, pName = fname, aLogical, aMother, false, 0);

  return aPhysVol;
}

G4bool G4PPVPlacement::IsMany() const
{
  return fmany; 
}

G4int G4PPVPlacement::GetCopyNo() const
{
  return fcopyNo;
}


void  G4PPVPlacement::SetCopyNo(G4int newCopyNo)
{
  fcopyNo= newCopyNo;
}

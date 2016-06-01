// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPhysicalVolume.cc,v 2.1 1998/07/12 02:55:41 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// class G4VPhysicalVolume Implementation

#include "G4VPhysicalVolume.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

// Constructor: init parameters and register in Store
G4VPhysicalVolume::G4VPhysicalVolume(G4RotationMatrix *pRot,
				     const G4ThreeVector &tlate,
				     const G4String& pName,
				     G4LogicalVolume* pLogical,
				     G4VPhysicalVolume* pMother)
{
    ftrans=tlate;
    frot=pRot;
    SetName(pName);
    SetLogicalVolume(pLogical);
    SetMother(pMother);
    if (pMother) pMother->GetLogicalVolume()->AddDaughter(this);
    G4PhysicalVolumeStore::Register(this);
}

// Destructor -  remove from Store
G4VPhysicalVolume::~G4VPhysicalVolume() 
{
    G4PhysicalVolumeStore::DeRegister(this);
}


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
// $Id: G4VPhysicalVolume.cc,v 1.4 2001-07-11 09:59:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4VPhysicalVolume Implementation
// 09.11.99 J.Apostolakis  Made GetObjectRotation() non-inline (it uses static)
//                         To be made obsolecent, in favour of more robust
//                           GetObjectRotationValue()

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

G4RotationMatrix* G4VPhysicalVolume::GetObjectRotation() const
{
  static G4RotationMatrix  aRotM; 
  static G4RotationMatrix  IdentityRM;  // Never changed (from "1")
  G4RotationMatrix* retval; 

  // Insure against frot being a null pointer
  if(frot)
    {
      aRotM= frot->inverse();
      retval= &aRotM;
    }
  else
    {
      retval= &IdentityRM;
    }
  return retval;
}

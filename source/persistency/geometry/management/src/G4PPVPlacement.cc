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
// $Id: G4PPVPlacement.cc,v 1.3 2001/07/11 10:02:20 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
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

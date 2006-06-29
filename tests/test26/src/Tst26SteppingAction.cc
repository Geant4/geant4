//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: Tst26SteppingAction.cc,v 1.6 2006-06-29 21:54:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26SteppingAction.hh"
#include "Tst26EventAction.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26SteppingAction::Tst26SteppingAction(Tst26EventAction* evt)
:Tst26Event(evt)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26SteppingAction::~Tst26SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.) return;

  const G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  const G4LogicalVolume* lv = pv->GetLogicalVolume();
  //  const G4MaterialCutsCouple* couple = aStep->GetTrack()->GetMaterialCutsCouple();
  //  G4int idx = couple->GetIndex();
  //  const G4Material* mat = couple->GetMaterial();
  //  G4cout << "Step in " << mat->GetName() << "   idx= " << idx << " edep= " << edep << G4endl;
  G4int volumeIndex = -1;
  
  G4int copyNo = pv->GetCopyNo();
  G4String name = lv->GetName();
  if(name == "Ecal") volumeIndex = 0;
  else if(name == "Abs1") volumeIndex = 1;
  else if(name == "Abs2") volumeIndex = 2;
  else if(name == "Abs3") volumeIndex = 3;
  else if(name == "Abs4") volumeIndex = 4;
  else if(name == "Vert") volumeIndex = 5;
  // G4cout << "     vIndx= " << volumeIndex << " copyNo= " << copyNo << G4endl;
  if(volumeIndex>=0) Tst26Event->AddEnergy(edep, volumeIndex, copyNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



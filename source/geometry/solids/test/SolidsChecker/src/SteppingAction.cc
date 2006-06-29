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
// $Id: SteppingAction.cc,v 1.3 2006-06-29 18:54:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm9: Crystal calorimeter
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

#include "SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"


/////////////////////////////////////////////////////////

SteppingAction::SteppingAction()
  //: theHisto(HistoManager::GetPointer())
{}

////////////////////////////////////////////////////////////

SteppingAction::~SteppingAction()
{}

/////////////////////////////////////////////////////////////////

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // theHisto->AddStep();
  // G4double edep = aStep->GetTotalEnergyDeposit();
  //  if(edep == 0.) return;

  const G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  const G4LogicalVolume* lv = pv->GetLogicalVolume();
  //  const G4MaterialCutsCouple* couple = aStep->GetTrack()->GetMaterialCutsCouple();
  //  G4int idx = couple->GetIndex();
  //  const G4Material* mat = couple->GetMaterial();
  //  G4cout << "Step in " << mat->GetName() << "   idx= " << idx << " edep= " << edep << G4endl;
  // G4int volumeIndex = -1;

  // G4int copyNo = pv->GetCopyNo();
  G4String name = lv->GetName();
  // G4cout<<name<<" "<<aStep->GetStepLength()<<"  ";
  if(name == "aVolume_L")
  {
    //  G4cout<<"l = "<<aStep->GetStepLength()<< G4endl;
  }
  // G4cout << "     vIndx= " << volumeIndex << " copyNo= " << copyNo << G4endl;
  // if(volumeIndex>=0) theHisto->AddEnergy(edep, volumeIndex, copyNo);
}

//////////////////////////////////////////////////////////////////



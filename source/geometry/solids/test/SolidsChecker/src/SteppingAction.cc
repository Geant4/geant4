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
// $Id: SteppingAction.cc,v 1.1 2004-11-26 16:58:10 grichine Exp $
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

  G4int copyNo = pv->GetCopyNo();
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



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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
:G4UserSteppingAction(),detector(det),eventAct(evt)
{
  first = true;
  lvol_world = lvol_module = lvol_layer = lvol_fiber = 0;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step )
{ 
 //some initialisation
 // 
 if (first) {
   lvol_world  = detector->GetLvolWorld();
   lvol_module = detector->GetLvolModule();   
   lvol_layer  = detector->GetLvolLayer();
   lvol_fiber  = detector->GetLvolFiber();      
   first = false;   
 }
 
 //if no edep, return
 //
 G4double edep = step->GetTotalEnergyDeposit();
 ///if (edep == 0.) return;
 
 //locate point in geometry
 //
 G4int iModule = 0;
 G4int iLayer  = 0;
 G4int iFiber  = 0;
   
 G4TouchableHandle touch1 = step->GetPreStepPoint()->GetTouchableHandle(); 
 G4LogicalVolume* lvol = touch1->GetVolume()->GetLogicalVolume();
  
      if (lvol == lvol_world) return;
 else if (lvol == lvol_module) { iModule = touch1->GetCopyNumber(0);}
 else if (lvol == lvol_layer)  { iLayer  = touch1->GetCopyNumber(0);
                                 iModule = touch1->GetCopyNumber(1);}
 else if (lvol == lvol_fiber)  { iFiber  = touch1->GetCopyNumber(0);
                                 iLayer  = touch1->GetCopyNumber(1);
                                 iModule = touch1->GetCopyNumber(2);}

 // sum edep
 //
 eventAct->SumDeStep(iModule, iLayer, iFiber, edep);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SteppingAction::BirksAttenuation(const G4Step* aStep)
{
 //Example of Birk attenuation law in organic scintillators.
 //adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
 //
 const G4Material* material = aStep->GetTrack()->GetMaterial();
 G4double birk1       = material->GetIonisation()->GetBirksConstant();
 G4double destep      = aStep->GetTotalEnergyDeposit();
 G4double stepl       = aStep->GetStepLength();  
 G4double charge      = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
 //
 G4double response = destep;
 if (birk1*destep*stepl*charge != 0.)
   {
     response = destep/(1. + birk1*destep/stepl);
   }
 return response;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


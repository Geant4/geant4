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
// $Id: SteppingAction.cc,v 1.9 2009-09-18 17:34:54 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4Step.hh"

#include "G4Geantino.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* run,
                               PrimaryGeneratorAction* prim, EventAction* evt,
			       HistoManager* hist)
:G4UserSteppingAction(),detector(det),runAct(run),primary(prim),eventAct(evt),
 histoManager(hist) 
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
   
   calorThickness  = detector->GetCalorThickness();
   calorSizeYZ     = detector->GetCalorSizeYZ();
   moduleThickness = detector->GetModuleThickness();   
   d1Pixel = detector->GetD1Pixel();
   d2Pixel = detector->GetD2Pixel();
   n1pxl   = detector->GetN1Pixels();
   n2pxl   = detector->GetN2Pixels();      
   n1shift = detector->GetN1Shift();   
      
   first = false;   
 }

 //locate point in geometry
 //  
 G4TouchableHandle touch1 = step->GetPreStepPoint()->GetTouchableHandle(); 
 G4LogicalVolume* lvol = touch1->GetVolume()->GetLogicalVolume();
 
 //if world, return
 //
 if (lvol == lvol_world) return;
 
 //sum nb of radiation length of calorimeter with geantino
 //
 G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();
 if (particle == G4Geantino::Geantino()) {
   G4double radl  = lvol->GetMaterial()->GetRadlen();
   G4double stepl = step->GetStepLength();
   eventAct->SumNbRadLength(stepl/radl);
 }
    
 //if no edep, return
 //
 G4double edep = step->GetTotalEnergyDeposit();
 if (edep == 0.) return;
 
 //locate position and compute pixel number
 //
 G4ThreeVector point1 = step->GetPreStepPoint()->GetPosition();
 G4int i1Module = (int) ((point1.x() + 0.5*calorThickness)/moduleThickness); 
 G4int i1Pixel  = (int) ((point1.x() + 0.5*calorThickness)/d1Pixel);
 if (i1Pixel < 0 ) i1Pixel = 0;
 if (i1Pixel >= n1pxl) i1Pixel = n1pxl - 1;
 G4double point1yz = point1.y();
 if (i1Module%2 != 0) point1yz = point1.z();
 G4int i2Pixel = (int) ((point1yz + 0.5*calorSizeYZ)/d2Pixel);
 if (i2Pixel < 0 ) i2Pixel = 0;
 if (i2Pixel >= n2pxl) i2Pixel = n2pxl - 1;
  
 G4int  iPixel = i1Pixel*n1shift + i2Pixel;
  
 // sum total energy deposit
 //
 eventAct->SumTotalEnergy(iPixel, edep);         

 //in fiber ?
 //  
 if (lvol == lvol_fiber) {
   eventAct->SumVisibleEnergy(iPixel, edep);
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SteppingAction::BirksAttenuation(const G4Step* aStep)
{
 //Example of Birk attenuation law in organic scintillators.
 //adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
 //
 G4Material* material = aStep->GetTrack()->GetMaterial();
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


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
// $Id: SteppingAction.cc,v 1.28 2008-05-29 16:59:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4Step.hh"
#include "G4Positron.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* run,
                               EventAction* evt, HistoManager* hist)
:G4UserSteppingAction(),detector(det),runAct(run),eventAct(evt),
 histoManager(hist) 
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //track informations
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();   
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition(); 
    
  //if World, return
  //
  G4VPhysicalVolume* volume = prePoint->GetTouchableHandle()->GetVolume();    
  //if sum of absorbers do not fill exactly a layer: check material, not volume.
  G4Material* mat = volume->GetLogicalVolume()->GetMaterial();
  if (mat == detector->GetWorldMaterial()) return; 
 
  //here we are in an absorber. Locate it
  //
  G4int absorNum  = prePoint->GetTouchableHandle()->GetCopyNumber(0);
  G4int layerNum  = prePoint->GetTouchableHandle()->GetCopyNumber(1);
         
  // collect energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  // collect step length of charged particles
  G4double stepl = 0.;
  if (particle->GetPDGCharge() != 0.) stepl = aStep->GetStepLength();
  
  //  G4cout << "Nabs= " << absorNum << "   edep(keV)= " << edep << G4endl;
  
  // sum up per event
  eventAct->SumEnergy(absorNum,edep,stepl);
  
  //longitudinal profile of edep per absorber
  if (edep>0.) histoManager->FillHisto(MaxAbsor+absorNum, 
				       G4double(layerNum+1), edep);
  
  //energy flow
  //
  // unique identificator of layer+absorber
  G4int Idnow = (detector->GetNbOfAbsor())*layerNum + absorNum;
  G4int plane;
  //
  //leaving the absorber ?
  if (endPoint->GetStepStatus() == fGeomBoundary) {
    G4ThreeVector position  = endPoint->GetPosition();
    G4ThreeVector direction = endPoint->GetMomentumDirection();
    G4double sizeYZ = 0.5*detector->GetCalorSizeYZ();       
    G4double Eflow = endPoint->GetKineticEnergy();
    if (particle == G4Positron::Positron()) Eflow += 2*electron_mass_c2;
    if ((std::abs(position.y()) >= sizeYZ) || (std::abs(position.z()) >= sizeYZ)) 
                                  runAct->sumLateralEleak(Idnow, Eflow);
    else if (direction.x() >= 0.) runAct->sumEnergyFlow(plane=Idnow+1, Eflow);
    else                          runAct->sumEnergyFlow(plane=Idnow,  -Eflow);    
  }   

////  example of Birk attenuation
///G4double destep   = aStep->GetTotalEnergyDeposit();
///G4double response = BirksAttenuation(aStep);
///G4cout << " Destep: " << destep/keV << " keV"
///       << " response after Birks: " << response/keV << " keV" << G4endl;
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


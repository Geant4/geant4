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
// $Id: SteppingAction.cc,v 1.20 2005/05/18 15:28:37 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"
#include "G4Positron.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* run,
                               EventAction* evt, HistoManager* hist)
:G4UserSteppingAction(),detector(det),runAct(run),eventAct(evt),
 histoManager(hist) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //initialize Energy flow
  if (aStep->GetTrack()->GetCurrentStepNumber() == 1) Idold = 0;
      
  //if World, return
  //
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4VPhysicalVolume* volume = prePoint->GetPhysicalVolume();    
  //if sum of absorbers do not fill exactly a layer: check material, not volume.
  G4Material* mat = volume->GetLogicalVolume()->GetMaterial();
  if (mat == detector->GetWorldMaterial()) return;
 
  //here we are in an absorber. Locate it
  //
  G4int absorNum  = volume->GetCopyNo();
  G4int layerNum  = prePoint->GetTouchable()->GetReplicaNumber(1);
  
  //track informations 
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4Track*     track    = aStep->GetTrack();
  const G4ParticleDefinition* particle = track->GetDefinition(); 
     
  // collect energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  // collect step length of charged particles
  G4double stepl = 0.;
  if (particle->GetPDGCharge() != 0.) stepl = aStep->GetStepLength();
    
  // sum up per event
  eventAct->SumEnergy(absorNum,edep,stepl);
  
  //longitudinal profile of edep per absorber
  if (edep>0.) histoManager->FillHisto(MaxAbsor+absorNum, layerNum+1., edep);
  
  //energy flow
  //
  // unique identificator of layer+absorber
  G4int Idnow = (detector->GetNbOfAbsor())*layerNum + absorNum;
  if (track->GetCurrentStepNumber() == 1) Idold = Idnow;     // track begins
  //
  //first step in the absorber ?
  G4int plane=0;
  if (prePoint->GetStepStatus() == fGeomBoundary) {
    G4double Eflow = prePoint->GetKineticEnergy();
    if (particle == G4Positron::Positron()) Eflow += 2*electron_mass_c2;  
    if      (Idnow > Idold) runAct->sumForwEflow(plane=Idnow,   Eflow);
    else if (Idnow < Idold) runAct->sumBackEflow(plane=Idnow+1, Eflow);
    Idold = Idnow;
  }   
  
  //last plane : forward leakage
  G4int Idlast = (detector->GetNbOfAbsor())*(detector->GetNbOfLayers()); 
  if ((Idnow == Idlast) && (endPoint->GetStepStatus() == fGeomBoundary)) {
    G4double Eflow = endPoint->GetKineticEnergy();
    if (particle == G4Positron::Positron()) Eflow += 2*electron_mass_c2;  
    runAct->sumForwEflow(plane=Idlast+1, Eflow);
  }

////  example of Birk attenuation
////  G4double destep   = aStep->GetTotalEnergyDeposit();
////  G4double response = BirkAttenuation(aStep);
////  G4cout << " Destep: " << destep/keV << " keV"
////         << " response after Birk: "  << response/keV << " keV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SteppingAction::BirkAttenuation(const G4Step* aStep)
{
 //Example of Birk attenuation law in organic scintillators.
 //adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
 //
 const G4String myMaterial = "Scintillator";
 const G4double birk1 = 0.013*g/(MeV*cm2);
 //
 G4double destep      = aStep->GetTotalEnergyDeposit();
 G4Material* material = aStep->GetTrack()->GetMaterial();
 G4double charge      = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
 //
 G4double response = destep;
 if ((material->GetName()==myMaterial)&&(charge!=0.))
   {
     G4double correction =
     birk1*destep/((material->GetDensity())*(aStep->GetStepLength()));
     response = destep/(1. + correction);
   }
 return response;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


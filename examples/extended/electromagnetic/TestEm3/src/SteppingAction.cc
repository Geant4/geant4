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
// $Id: SteppingAction.cc,v 1.14 2004/11/23 14:05:31 maire Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"
#include "G4Positron.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt,
                               HistoManager* hist)
:G4UserSteppingAction(),detector(det),eventAct(evt),histoManager(hist) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* endPoint = aStep->GetPostStepPoint();
  
  //if World, returns
  //
  G4VPhysicalVolume* volume = prePoint->GetPhysicalVolume();
  //if sum of absorbers do not fill exactly a layer: check material, not volume.
  G4Material* mat = volume->GetLogicalVolume()->GetMaterial();
  if (mat == detector->GetWorldMaterial()) return;
  
  //locate the absorber
  //
  G4int absorNum  = volume->GetCopyNo();
  G4int layerNum  = prePoint->GetTouchable()->GetReplicaNumber(1);
    
  ////G4int absorNum  = prePoint->GetTouchable()->GetCopyNumber( );
  ////G4int layerNum  = prePoint->GetTouchable()->GetCopyNumber(1);
   
  // collect energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  // collect step length of charged particles
  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();
    
  // sum up per event
  eventAct->SumEnergy(absorNum,edep,stepl);
  
  //longitudinal profile of edep per absorber
  histoManager->FillHisto(MaxAbsor+absorNum, layerNum+1., edep);
  
  //energy flow
  //
  //leaving an absorber ?  in forward direction ?
  const G4Track* track = aStep->GetTrack();
  if ((endPoint->GetProcessDefinedStep()->GetProcessName() == "Transportation")
      && (track->GetMomentumDirection().x() > 0.)) {
      G4int planNum = 1 + (detector->GetNbOfAbsor())*layerNum + absorNum;
      G4double EnLeaving = track->GetKineticEnergy();
      if (track->GetDefinition() == G4Positron::Positron())
       EnLeaving += 2*electron_mass_c2;
      G4int ih = 2*MaxAbsor + 1;
      if (track->GetTrackID() != 1) ih += 1;
      histoManager->FillHisto(ih, (G4double)planNum, EnLeaving);
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


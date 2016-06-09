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
// $Id: SteppingAction.cc,v 1.9 2004/01/21 17:29:27 maire Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Track.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
:G4UserSteppingAction(),detector(det),eventAct(evt) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  const G4Track* track = aStep->GetTrack();
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4Material* mat = volume->GetLogicalVolume()->GetMaterial();
  if (mat == detector->GetWorldMaterial()) return;
  
  G4int absorNo = volume->GetCopyNo();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();

  G4double stepl = 0.;
  if (track->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();
    
  // sum up per event
  eventAct->SumEnergy(absorNo,edep,stepl);
  
  // energy leaving each volume
  if (track->GetNextVolume() != volume) {
    G4double EnLeaving = track->GetKineticEnergy();
    if (track->GetDefinition() == G4Positron::Positron())
       EnLeaving += 2*electron_mass_c2;
    eventAct->SumEnLeaving(absorNo,EnLeaving);
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


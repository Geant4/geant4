// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3SteppingAction.cc,v 1.3 2001-02-20 12:34:44 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em3SteppingAction.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3SteppingAction::Em3SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3SteppingAction::~Em3SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3SteppingAction::UserSteppingAction(const G4Step* aStep)
{
////  G4double destep   = aStep->GetTotalEnergyDeposit();
////  G4double response = BirkAttenuation(aStep);
////  G4cout << " Destep: " << destep/keV << " keV"   
////         << " response after Birk: "  << response/keV << " keV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double Em3SteppingAction::BirkAttenuation(const G4Step* aStep)
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ios.hh"

#include "ZIIISteppingAction.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh" 
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"

#include "g4std/vector"

extern G4bool drawEvent;

extern G4std::vector<G4String> Particles;
extern G4std::vector<G4double> Energies;
extern G4std::vector<G4double> Weights;
extern G4std::vector<G4double> Times;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ZIIISteppingAction::ZIIISteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ZIIISteppingAction::~ZIIISteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ZIIISteppingAction::UserSteppingAction(const G4Step* fStep)
{ 
  const G4SteppingManager* pSM = fpSteppingManager;
  G4Track* fTrack = pSM->GetTrack();
  //G4Step* fStep = pSM->GetStep();
  //G4int TrackID = fTrack->GetTrackID();
  G4int StepNo = fTrack->GetCurrentStepNumber();
  if(StepNo >= 10000) fTrack->SetTrackStatus(fStopAndKill);
  
  //  cout << fTrack->GetGlobalTime()/s <<"  "<< fTrack->GetLocalTime()/s << " " <<fTrack->GetProperTime()/s << endl;
  //cout << fStep->GetPreStepPoint()->GetGlobalTime() /s << endl; 

  if (StepNo == 1) {
    Particles.push_back ( fTrack->GetDefinition()->GetParticleName() );
    Energies.push_back ( fStep->GetPreStepPoint()->GetKineticEnergy()/keV );
    Weights.push_back ( fStep->GetPreStepPoint()->GetWeight() );
    Times.push_back((fStep->GetPreStepPoint()->GetGlobalTime() - fStep->GetPreStepPoint()->GetLocalTime()) / s );
    drawEvent = true;
  }

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if (pVVisManager) {

        //----- Get the Stepping Manager
        const G4SteppingManager* pSM = fpSteppingManager;

        //----- Define a line segment 
        G4Polyline polyline;

        G4String name = pSM->GetTrack()->GetDefinition()->GetParticleName();
//        G4double charge = pSM->GetTrack()->GetDefinition()->GetPDGCharge();
//        G4double mass = pSM->GetTrack()->GetDefinition()->GetPDGMass();
        G4Colour colour;
        if  (name == "neutron") colour = G4Colour(1., 0., 0.);
        else if (name == "gamma" ) colour = G4Colour(0., 0., 1.);
        else                  colour = G4Colour(0., 1., 0.);
        G4VisAttributes attribs(colour);
        polyline.SetVisAttributes(attribs);
        G4Point3D * start = new G4Point3D(pSM->GetStep()->GetPreStepPoint()->GetPosition() );
        G4Point3D * end = new G4Point3D(pSM->GetStep()->GetPostStepPoint()->GetPosition() );
        polyline.push_back(*start);
        polyline.push_back(*end);
        delete start;
        delete end;

        //----- Call a drawing method for G4Polyline 
        pVVisManager -> Draw(polyline); // Step 2 of the Clear-Draw-Show process
      





  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




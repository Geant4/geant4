//  XrayTelSteppingAction.cc

#include "G4ios.hh"
#include <assert.h>
#include <fstream.h>
#include <iomanip.h>
#include <iostream.h>

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "XrayTelSteppingAction.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"

#include "g4std/vector"

extern G4bool drawEvent;
extern G4std::vector<G4String> EnteringParticles;
extern G4std::vector<G4double> EnteringEnergy;
extern G4std::vector<G4ThreeVector> EnteringDirection;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelSteppingAction::XrayTelSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelSteppingAction::~XrayTelSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelSteppingAction::UserSteppingAction(const G4Step*)
{
  const G4SteppingManager* pSM = fpSteppingManager;
  G4Track* fTrack = pSM->GetTrack();
  G4Step* fStep = pSM->GetStep();
  G4int TrackID = fTrack->GetTrackID();
  G4int StepNo = fTrack->GetCurrentStepNumber();


  if(StepNo >= 10000) fTrack->SetTrackStatus(fStopAndKill);

  G4String volName; 
  if ( fTrack->GetVolume() ) 
    volName =  fTrack->GetVolume()->GetName(); 
  G4String nextVolName;
  if ( fTrack->GetNextVolume() ) 
    nextVolName =  fTrack->GetNextVolume()->GetName();
 
  Hep3Vector pos = fTrack->GetPosition();

  //--- Entering Detector
  if(volName != "Detector_P" && nextVolName == "Detector_P") {
    EnteringParticles.push_back ( fTrack->GetDefinition()->GetParticleName() );
    EnteringEnergy.push_back ( fTrack->GetKineticEnergy() );
    EnteringDirection.push_back (pos);

    drawEvent = true;
  }
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




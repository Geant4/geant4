// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelSteppingAction.cc                        *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope Physics list
// - Based on Chandra and XMM models
// 
//
// **********************************************************************

#include "G4ios.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"

#include "globals.hh"

#include <assert.h>
#include "g4std/fstream"
#include "g4std/iomanip"
#include "g4std/iostream"
#include "g4std/vector"

#include "XrayTelSteppingAction.hh"
//#include "XrayTelHistogram.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelSteppingAction::XrayTelSteppingAction(
					     G4std::vector<G4double*>* enEnergy, 
					     G4std::vector<G4ThreeVector*>* enDirect,
					     G4bool* dEvent)
  : EnteringEnergy(enEnergy),
    EnteringDirection(enDirect),drawEvent(dEvent)
{;}

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
 
  G4ThreeVector pos = fTrack->GetPosition();

  //--- Entering Detector
  if(volName != "Detector_P" && nextVolName == "Detector_P") {

    EnteringEnergy->push_back ( new G4double (fTrack->GetKineticEnergy()) );
    EnteringDirection->push_back (new G4ThreeVector (pos));

    // now we want to do some analysis at this step ... 
    // call back to the analysis-manger to do the analysis ...
    //    histoManager->analyze(pSM);

    *drawEvent = true;
  }
}








//  XrayTelSteppingAction.cc

#include "G4ios.hh"
#include <assert.h>
#include <fstream.h>
#include <iomanip.h>
#include <iostream.h>

//#include "CLHEP/Hist/TupleManager.h"
//#include "CLHEP/Hist/HBookFile.h"
//#include "CLHEP/Hist/Histogram.h"
//#include "CLHEP/Hist/Tuple.h"

#include "XrayTelHistogram.hh"

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

XrayTelSteppingAction::XrayTelSteppingAction(XrayTelHistogram* hMgr)
 : histoManager(hMgr)
{
}

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

    IHistogram1D * kinE;
    if (histoManager == 0) {
      G4cerr << "histomanger not there !" << G4endl;
    } else {
      vector<IHistogram1D *> * hlist = histoManager->getH1DList();
      if (hlist == 0) {
	G4cerr << " 1d list not found " << G4endl;
      } else {
	kinE = (*hlist)[0]; // 1D-histo # 0 is kinetic engergy
       	if (kinE == 0) {
	  G4cerr << "XrayTelSteppingAction::UserSteppingAction> could not find histo for kinetic energy !" << G4endl;
	} else {
	  G4cerr << "filling histo" << G4endl;
	  kinE->fill(fTrack->GetKineticEnergy());
	}
      }
    }
    
    IHistogram2D * posHist;
    posHist = (*(histoManager->getH2DList()))[0];// 2D-histo # 0 is pos
    if (posHist == 0) {
      G4cerr << "XrayTelSteppingAction::UserSteppingAction> could not find histo for position !" << G4endl;
    } else {
      posHist->fill(pos.y(), pos.z());
      histoManager->plot(posHist);
      histoManager->getPlotter()->refresh();
      histoManager->getPlotter()->psPrint("posplot.ps");
    }
    
    drawEvent = true;
  }
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




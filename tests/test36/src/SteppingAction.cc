#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4Track.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* arun)
:G4UserSteppingAction(),detector(det),runaction(arun) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  const G4Track* track = aStep->GetTrack();
  G4VPhysicalVolume* volume = track->GetVolume();
  G4double edep = aStep->GetTotalEnergyDeposit();

  if ((edep!=0.)&&(volume -> GetName() == "LayerMedium 1"))
  {// sum per layer
    runaction->SumEnergy1(volume->GetCopyNo(),edep);
  // sum per medium  
    runaction -> AddEnergy1(edep);
  }

  if ((edep!=0.)&&(volume -> GetName() == "LayerMedium 2"))
  { runaction->SumEnergy2(volume->GetCopyNo(),edep);
    runaction -> AddEnergy2(edep);
  }
  
  if ((edep!=0.)&&(volume -> GetName() == "LayerMedium 3"))
  { runaction->SumEnergy3(volume->GetCopyNo(),edep);
    runaction -> AddEnergy3(edep);
  }
  
  
}


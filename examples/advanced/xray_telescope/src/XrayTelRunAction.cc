//  XrayTelRunAction.cc

#include "XrayTelRunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include <fstream.h>

#include "g4std/vector"

extern G4bool drawEvent;
extern G4std::vector<G4String> EnteringParticles;
extern G4std::vector<G4double> EnteringEnergy;
extern G4std::vector<G4ThreeVector> EnteringDirection;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelRunAction::XrayTelRunAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelRunAction::~XrayTelRunAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int RunN = aRun->GetRunID();
  if ( RunN % 1000 == 0 ) 
    G4cout << "### Run : " << RunN << endl;

  if (G4VVisManager::GetConcreteInstance()) {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/clear/view");
      UI->ApplyCommand("/vis/draw/current");
  } 
  
  EnteringParticles.clear();
  EnteringEnergy.clear();
  EnteringDirection.clear();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelRunAction::EndOfRunAction(const G4Run* )
{
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/show/view");

  ofstream outscat("detector.hist", ios::app);

  cout << "End of Run summary" << endl << endl;

  G4double TotEnteringEnergy = 0.0;
  for (int i=0;i<EnteringParticles.size();i++)
	TotEnteringEnergy +=  EnteringEnergy[i];
  cout << "Total Entering Detector : " << EnteringParticles.size()  << endl;
  cout << "Total Entering Detector Energy : " << TotEnteringEnergy  << endl;

  for (i=0;i<EnteringParticles.size();i++) {
    outscat << "  "
	    << EnteringEnergy[i]
            << "  "
            << EnteringDirection[i].x()
	    << "  "
            << EnteringDirection[i].y()
	    << "  "
            << EnteringDirection[i].z()
	    << endl;
  }
  outscat.close();
}



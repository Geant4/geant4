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
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "ZIIIRunActionMessenger.hh"
#include "ZIIIRunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "g4std/fstream"
#include "g4std/iomanip"
#include "g4std/vector"

extern G4std::vector<G4String> Particles;
extern G4std::vector<G4double> Energies;
extern G4std::vector<G4double> Weights;
extern G4std::vector<G4double> Times;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ZIIIRunAction::ZIIIRunAction()
  : fileName("rdmex2.log")
{
  runMessenger = new ZIIIRunActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ZIIIRunAction::~ZIIIRunAction()
{
  delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ZIIIRunAction::BeginOfRunAction(const G4Run* aRun)
{
 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 

  Particles.clear();
  Energies.clear();
  Weights.clear();
  Times.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ZIIIRunAction::EndOfRunAction(const G4Run* aRun)
{
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  G4std::ofstream outscat(fileName, G4std::ios::app);

  for (G4int i=0; i<Particles.size();i++) {

    outscat 
      << G4std::setiosflags(G4std::ios::fixed)
      << G4std::setprecision(5)
      << G4std::setiosflags(G4std::ios::right)
      << G4std::setw(12)
      << Energies[i]
      << G4std::setw(12)<<G4std::setprecision(5) 
      << G4std::setiosflags(G4std::ios::scientific)
      << G4std::setiosflags(G4std::ios::right)
      << Weights[i]
      << G4std::setw(12)<<G4std::setprecision(5)
      << G4std::setiosflags(G4std::ios::scientific)
      << G4std::setiosflags(G4std::ios::right)
      << Times[i] << "     "
      << Particles[i]
      << G4endl ;    
  }
  outscat << G4endl;
  outscat.close();



  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

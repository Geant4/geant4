// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst18RunAction.cc,v 1.6 2000-06-15 09:52:40 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst18RunActionMessenger.hh"
#include "Tst18RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
//#include "Histo.hh"

#include "g4std/fstream"
#include "g4std/iomanip"
#include "g4std/vector"

//using namespace std;

//extern G4String filename;

extern G4std::vector<G4String> Particles;
extern G4std::vector<G4double> Energies;
extern G4std::vector<G4double> Weights;
extern G4std::vector<G4double> Times;
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst18RunAction::Tst18RunAction()
  : fileName("rdmex2.log")
{
  runMessenger = new Tst18RunActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst18RunAction::~Tst18RunAction()
{
  delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst18RunAction::BeginOfRunAction(const G4Run* aRun)
{
 
  G4int RunN = aRun->GetRunID();
  if ( RunN % 1000 == 0 ) 
    G4cout << "### Run : " << RunN << G4endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/clear/view");
      UI->ApplyCommand("/vis/draw/current");
    } 
  
  Particles.clear();
  Energies.clear();
  Weights.clear();
  Times.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst18RunAction::EndOfRunAction(const G4Run* )
{
  
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/show/view");
  //  ofstream outscat(fileName, ios::app);

  for (G4int i=0; i<Particles.size();i++) {

    //    outscat 
    G4cout
      << G4std::setiosflags(G4std::ios::fixed)
      << G4std::setprecision(3)
      << G4std::setiosflags(G4std::ios::right)
      << G4std::setw(12)
      << Energies[i]
      << G4std::setw(12)<<G4std::setprecision(4) 
      << G4std::setiosflags(G4std::ios::scientific)
      << G4std::setiosflags(G4std::ios::right)
      << Weights[i]
      << G4std::setw(12)<<G4std::setprecision(4)
      << G4std::setiosflags(G4std::ios::scientific)
      << G4std::setiosflags(G4std::ios::right)
      << Times[i] << "     "
      << Particles[i]
      << G4endl ;    
  }
  G4cout<< G4endl;
  //  outscat << G4endl;
  //outscat.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





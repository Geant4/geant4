#include "G4AtomicTransitionManager.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"

int main() {


  G4int Z;
  G4int subShell;

  G4cout << "Enter Z " << G4endl;
  G4cin >> Z;

  G4AtomicTransitionManager* transManager = G4AtomicTransitionManager::Instance();

  G4cout << "Number of possible shells: " << transManager->NumberOfShells(Z);
  G4cout << "Select SubShell: " << G4endl;
  G4cin >> subShell;

  G4cout << "Testing G4AtomicShell" << G4endl;
  G4cout << "Primary Shell: " << transManager->Shell(Z,subShell)->ShellId() << G4endl;
  G4cout << "SubShell binding energy: " << transManager->Shell(Z,subShell)->BindingEnergy() << G4endl;
  
  G4std::vector<G4double> transEnergies = transManager->Shell(Z,subShell)->TransitionEnergies();
  G4std::vector<G4int> transIds = transManager->Shell(Z,subShell)->TransSubShellIdentifiers();
  G4std::vector<G4double> transProbs = transManager->Shell(Z,subShell)->TransitionProbabilities();

  for (G4int trans=0; trans<transIds.size(); trans++) {

    G4cout << "SubShell transition Start from TransSubShellIdentifiers: " << transIds[trans] << G4endl;
    G4cout << "SubShell transition Start from TransitionIdentifier: " << 
      transManager->Shell(Z,subShell)->TransitionIdentifier(trans) << G4endl;

    G4cout << "SubShell transition Start from TransitionEnergies: " << transEnergies[trans] << G4endl;
    G4cout << "SubShell transition Start from TransitionEnergy: " << 
      transManager->Shell(Z,subShell)->TransitionEnergy(trans) << G4endl;

    G4cout << "SubShell transition Start from TransitionProbabilities: " << transProbs[trans] << G4endl;
    G4cout << "SubShell transition Start from TransitionProbability: " << 
      transManager->Shell(Z,subShell)->TransitionProbability(trans) << G4endl;

  }

  G4cout << "Testing G4AtomicTransitionManager" << G4endl;

  G4cout << "Total number of SubShells: " << transManager->NumberOfShells(Z) <<G4endl;

  G4cout << "Total Probability that a radiative transition occurs: " <<
             transManager->TotalRadiativeTransitionProbability(Z,subShell) <<G4endl;

  G4cout << "Total Probability that a NON radiative transition occurs: " <<
             transManager->TotalNonRadiativeTransitionProbability(Z,subShell) <<G4endl;

}

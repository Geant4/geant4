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

  G4cout << "Number of shells: " << transManager->NumberOfShells(Z)<<G4endl;



  G4cout<< "Number of Fluo vacancies: "<<transManager->NumberOfReachableShells(Z)<<G4endl;



  
  G4cout << "Testing G4AtomicShell" << G4endl;
  G4cout << "Select the index of the subshell whose binding energy you need: " << G4endl;
  G4cin >> subShell;
  G4cout << "Primary Shell: " << transManager->Shell(Z,subShell)->ShellId() << G4endl;
  G4cout << "SubShell binding energy: " << transManager->Shell(Z,subShell)->BindingEnergy() << G4endl;
  
  G4cout << "Testing G4AtomicShell" << G4endl;
  G4int shellIndex;


  G4cout <<" Select the index of the Fluo vacancy "<<G4endl;
  G4cin>> shellIndex; 
  
  G4cout << "Testing G4FluoTransition "<<G4endl;
  G4std::vector<G4double> transEnergies = transManager->ReachableShell(Z,shellIndex)->TransitionEnergies();
  G4std::vector<G4int> transIds = transManager->ReachableShell(Z,shellIndex)->OriginatingShellIds();
  G4std::vector<G4double> transProbs = transManager->ReachableShell(Z,shellIndex)->TransitionProbabilities();

  for (G4int trans=0; trans<transIds.size(); trans++) {

    G4cout << "The transition starts from the shell: " << transIds[trans] << G4endl;
    G4cout << "The transition starts from the shell: " << 


      transManager->ReachableShell(Z,shellIndex)->OriginatingShellId(trans) << G4endl;


    G4cout << " Transition energy: " << transEnergies[trans] << G4endl;
    G4cout << "Transition energy: " << 
      transManager->ReachableShell(Z,shellIndex)->TransitionEnergy(trans) << G4endl;

    G4cout << "Transition probability: " << transProbs[trans] << G4endl;
    G4cout << "Transition probability: " << 
      transManager->ReachableShell(Z,shellIndex)->TransitionProbability(trans) << G4endl;

  }

  /* ===Atention=== the data given out by the following code usually doesn't
     corespond to the data of the preceding one, I.E. the vacancy used for the testing of auger
     effect isn't the same of that used for Fluorescence.
  */



  G4cout<< "Number of Auger vacancies: "<<transManager->NumberOfReachableAugerShells(Z)<<G4endl;
  G4cout <<" Select the index of the Auger vacancy "<<G4endl;
  G4int augerVacancyIndex(0);
  G4cin>> augerVacancyIndex; 


  G4cout << "Testing G4AugerTransition "<<G4endl;

  const G4AugerTransition* augerTransition = transManager->ReachableAugerShell(Z, augerVacancyIndex);
  
  const G4std::vector<G4int> augerTransIds = *(augerTransition->TransitionOriginatingShellIds());
  for (G4int transIndex = 0; transIndex<(augerTransIds.size() ); transIndex++) {
 
    G4cout << "Questo e' augerTransIds[transIndex]: "<< augerTransIds[transIndex] << G4endl;

    G4std::vector<G4double> augerProbs = *augerTransition->AugerTransitionProbabilities(augerTransIds[transIndex]);
    G4std::vector<G4int> augerIds = *augerTransition->AugerOriginatingShellIds(augerTransIds[transIndex]);
    G4std::vector<G4double> augerEnergies = *augerTransition->AugerTransitionEnergies(augerTransIds[transIndex]);
    
    // we r comparing the vectors elements given by G4AugerTransition with the values of the single functions.
    
    G4cout << "The transition starts from the shell: " << augerTransIds[transIndex] << G4endl;

    G4int transId = augerTransition->TransitionOriginatingShellId(transIndex);

    G4cout << "The transition starts from the shell: " << transId << G4endl;
    
    for (G4int trans=0; trans<augerIds.size(); trans++) {
      
      G4cout << "The auger electron starts from the shell: " << augerIds[trans] << G4endl;
      G4cout << "The auger electron starts from the shell: " << 


	augerTransition->AugerOriginatingShellId(trans, transId) << G4endl;


      
      G4cout << " Transition energy: " << augerEnergies[trans] << G4endl;
      G4cout << "Transition energy: " << 
	augerTransition->AugerTransitionEnergy(trans, transId) << G4endl;
      
      G4cout << "Transition probability: " << augerProbs[trans] << G4endl;
      G4cout << "Transition probability: " << 
	 augerTransition->AugerTransitionProbability(trans, transId) << G4endl;
      
    }

  }


  G4cout << "Testing G4AtomicTransitionManager" << G4endl;

  G4cout << "Total number of SubShells: " << transManager->NumberOfShells(Z) <<G4endl;

  G4cout <<"Number of available shells: "<< transManager->NumberOfReachableShells(Z);

  G4cout << "Total Probability that a radiative transition occurs: " <<
             transManager->TotalRadiativeTransitionProbability(Z,shellIndex) <<G4endl;

  G4cout << "Total Probability that a NON radiative transition occurs: " <<
             transManager->TotalNonRadiativeTransitionProbability(Z,shellIndex) <<G4endl;

  G4cout << "END OF THE MAIN PROGRAM "<<G4endl;
}


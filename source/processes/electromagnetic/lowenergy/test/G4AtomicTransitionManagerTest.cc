//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "G4AtomicTransitionManager.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <vector>
#include <assert.h>

int main() {



  G4int i; 
  G4int zMin;
  G4int zMax;

  G4cout << "Enter Z " << G4endl;
  G4cin >> i;

  if (i == 0) 
    { 
    zMin = 6; 
    zMax = 100;
    }
  else if (i>0) 
    {
    zMin = i; 
    zMax = i;
    }  

  for (G4int Z=zMin; Z<=zMax; Z++) {

    G4AtomicTransitionManager* transManager = G4AtomicTransitionManager::Instance();

    G4int shellNumber = transManager->NumberOfShells(Z);
    G4cout << "Number of shells: " << shellNumber<<G4endl;
    

    G4cout << "Testing G4AtomicShell" << G4endl;
    G4int subShell ;
    G4int subShellIndexMax;
    G4int subShellIndexMin;
    if (i == 0) 
    { 
      subShellIndexMin = 0;
      subShellIndexMax = shellNumber;
    }
    else if (i > 0)
      {
    G4cout << "Select the index of the subshell whose binding energy you need: " << G4endl;
    G4cin >> subShell;
    subShellIndexMin = subShell;
    subShellIndexMax = subShell;
      }

    for (subShell= subShellIndexMin; subShell <= subShellIndexMax; subShell++) {

      G4cout << "Shell ID: " << transManager->Shell(Z,subShell)->ShellId() << G4endl;
      G4cout << "Shell binding energy: " << transManager->Shell(Z,subShell)->BindingEnergy() << G4endl;
    }

    G4int fluoVacNumber = transManager->NumberOfReachableShells(Z);
    
    G4cout<< "Number of Fluo vacancies: "<<fluoVacNumber<<G4endl;

    G4int shellIndex;
    G4int shellIndexMax;
    G4int shellIndexMin;
    if (i == 0) 
    { 
      shellIndexMin = 0;
      shellIndexMax = fluoVacNumber;
    }
    if (i >0)
      {
	G4cout <<" Select the index of the Fluo vacancy "<<G4endl;
	G4cin>> shellIndex; 
	shellIndexMin = shellIndex;
	shellIndexMax = shellIndex;     
      }

    for (shellIndex=shellIndexMin; shellIndex<=shellIndexMax; shellIndex++) {
      
      std::vector<G4double> transEnergies;
      
      const G4FluoTransition* transIt = transManager->ReachableShell(Z,shellIndex);
      
      G4cout << G4endl << "Testing G4FluoTransition "<<G4endl;
      if ( transIt ) {transEnergies = transManager->ReachableShell(Z,shellIndex)->TransitionEnergies();
	std::vector<G4int> transIds = transManager->ReachableShell(Z,shellIndex)->OriginatingShellIds();
	std::vector<G4double> transProbs = transManager->ReachableShell(Z,shellIndex)->TransitionProbabilities();
	
	
	for (G4int trans=0; trans<transIds.size(); trans++) {
	  
	  G4cout << G4endl << "Vacancy filled by e- from shell [vector]: " << transIds[trans] << G4endl;
	  G4cout << "Vacancy filled by e- from shell [single eval]: " << 
	    transManager->ReachableShell(Z,shellIndex)->OriginatingShellId(trans) << G4endl;
	  
	  G4cout << "Photon energy [vector]: " << transEnergies[trans] << G4endl;
	  G4cout << "Photon energy [single eval]: " << 
	    transManager->ReachableShell(Z,shellIndex)->TransitionEnergy(trans) << G4endl;
	  
	  G4cout << "Transition probability [vector]: " << transProbs[trans] << G4endl;
	  G4cout << "Transition probability [single eval]: " << 
	    transManager->ReachableShell(Z,shellIndex)->TransitionProbability(trans) << G4endl;
	  
	}
      }
    }
    /* ===Atention=== the data given out by the following code usually doesn't
       corespond to the data of the preceding one, I.E. the vacancy of index IDX for Auger 
       is isn't the same for Fluorescence.
    */


    G4int augerVacanciesNumber = transManager->NumberOfReachableAugerShells(Z);
    G4cout<< G4endl << "Number of Auger vacancies: "<< augerVacanciesNumber << G4endl;

    G4int augerVacancyIndex;
    G4int augerVacancyIndexMax;
    G4int augerVacancyIndexMin;
    if (i == 0) 
    { 
      augerVacancyIndexMin = 0;
      augerVacancyIndexMax = augerVacanciesNumber;
    }
    if (i >0)
      {

	G4cout <<" Select the index of the Auger vacancy "<<G4endl;
	G4int augerVacancyIndex(0);
	G4cin>> augerVacancyIndex; 
	augerVacancyIndexMin = augerVacancyIndex;
	augerVacancyIndexMax = augerVacancyIndex;     
      }
    
    for (augerVacancyIndex=augerVacancyIndexMin; augerVacancyIndex<=augerVacancyIndexMax; augerVacancyIndex++) {
      
      
      
      G4cout << G4endl << "Testing G4AugerTransition "<<G4endl;
      
      const G4AugerTransition* augerTransition = transManager->ReachableAugerShell(Z, augerVacancyIndex);
      
      const std::vector<G4int> augerTransIds = *(augerTransition->TransitionOriginatingShellIds());
      for (G4int transIndex = 0; transIndex<(augerTransIds.size() ); transIndex++) {
	
	G4cout << G4endl << "augerTransIds[transIndex]: "<< augerTransIds[transIndex] << G4endl;
	
	std::vector<G4double> augerProbs = *augerTransition->AugerTransitionProbabilities(augerTransIds[transIndex]);
	std::vector<G4int> augerIds = *augerTransition->AugerOriginatingShellIds(augerTransIds[transIndex]);
	std::vector<G4double> augerEnergies = *augerTransition->AugerTransitionEnergies(augerTransIds[transIndex]);
	
	// we r comparing the vectors elements given by G4AugerTransition with the values of the single functions.
	
	G4cout << "Vacancy filled by e- from the shell [vector]: " << augerTransIds[transIndex] << G4endl;
	
	G4int transId = augerTransition->TransitionOriginatingShellId(transIndex);
	
	G4cout << "Vacancy filled by e- from the shell [single eval]: " << transId << G4endl;
	
	for (G4int trans=0; trans<augerIds.size(); trans++) {
	  
	  G4cout << "The Auger e- starts from the shell [vector]: " << augerIds[trans] << G4endl;
	  G4cout << "The Auger e- starts from the shell [single eval]: " << 
	    
	    
	    augerTransition->AugerOriginatingShellId(trans, transId) << G4endl;
	  	  
	  G4cout << "Auger energy [vector]: " << augerEnergies[trans] << G4endl;
	  G4cout << "Auger energy [single eval]: " << 
	    augerTransition->AugerTransitionEnergy(trans, transId) << G4endl;
	  
	  G4cout << "Transition probability [vector]: " << augerProbs[trans] << G4endl;
	  G4cout << "Transition probability [single eval]: " << 
	    augerTransition->AugerTransitionProbability(trans, transId) << G4endl;
	  
	}
	
      }
    }
    
    G4cout << G4endl << "Testing G4AtomicTransitionManager" << G4endl;
    
    G4cout << "Total number of SubShells: " << transManager->NumberOfShells(Z) <<G4endl;
    
    G4cout <<"Number of available shells: "<< transManager->NumberOfReachableShells(Z)<< G4endl;
    


    for (shellIndex=shellIndexMin; shellIndex<=shellIndexMax; shellIndex++){

    G4cout << "Total Probability that a radiative transition occurs for shell "<<transManager->Shell(Z,shellIndex)->ShellId() << ": " << transManager->TotalRadiativeTransitionProbability(Z,shellIndex) <<G4endl;

  }
    for (shellIndex=augerVacancyIndexMin; shellIndex<=augerVacancyIndexMax; shellIndex++) {
          
      G4cout << "Total Probability that a NON radiative transition occurs for shell " << transManager->Shell(Z,shellIndex)->ShellId()<< ": " <<
      transManager->TotalNonRadiativeTransitionProbability(Z,shellIndex) <<G4endl;
    }
  }
  G4cout << G4endl << "END OF THE MAIN PROGRAM "<<G4endl;
}


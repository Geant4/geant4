// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4AtomicTransitionManager.cc
//
//      Author:        Alfonso Mantero (Alf@mailandnews.com)
// 
//      Creation date: 4 May 2001
// -------------------------------------------------------------------

#include "G4AtomicTransitionManager.hh"
#include "G4LowEnergyUtilities.hh"

//----------..................-------------------....................

// constructor

G4AtomicTransitionManager::G4AtomicTransitionManager() {

  BuildTable();

}

//----------..................-------------------....................
// this is destructor, private, because from outside i can only see instance()

G4AtomicTransitionManager::~G4AtomicTransitionManager() {

  if (instance) {

    G4std::map <G4int, G4std::vector<G4AtomicShell> >::iterator cur;
    for(cur = shellTable.begin(); cur != shellTable.end(); ++cur) {

      shellTable.erase(cur);
    }

    instance = 0;

  }
  

}

//....................000000000000000...................000000000000............00000000000000

G4int G4AtomicTransitionManager::NumberOfShells(G4int Z) {

  return shellTable[Z].size(); 

}

//....................000000000000000...................000000000000............00000000000000

// This returns the probability a fluorescence photon is emitted. It is the sum of the probabilities
// of emission for all possible transitions of the selected subshell

G4double G4AtomicTransitionManager::TotalRadiativeTransitionProbability(G4int Z, G4int subShellId) {
  
  G4double theRadTransProb(0);
 
  G4std::vector<G4int> subShellVect = Shell(Z,subShellId)->TransSubShellIdentifiers();

  for(G4int transId =0; transId != subShellVect.size(); transId++) {

    theRadTransProb += Shell(Z,subShellId)->TransitionProbability(transId);
  }

  return theRadTransProb;

}

//....................000000000000000...................000000000000............00000000000000

// since no auger effect is simulated yet, and first of all because probabilities r normalized 
// to one on all the transitions (radiative and non-radiative), total non radiative transition probability
// is computed as 1-(total radiative tarnsition probability)

G4double G4AtomicTransitionManager::TotalNonRadiativeTransitionProbability(G4int Z, G4int subShellId) {

  G4double nonRadProb = (1.0)-(TotalRadiativeTransitionProbability(Z, subShellId));
  return nonRadProb;

}
//....................000000000000000...................000000000000............00000000000000


G4AtomicShell* G4AtomicTransitionManager::Shell(G4int Z, G4int subShellId) {

  G4AtomicShell* tmp;
  *tmp = shellTable[Z][subShellId];

  return tmp;

}

//....................000000000000000...................000000000000............00000000000000

// this is a trick to get one only object of this type to exist at the same time

G4AtomicTransitionManager* G4AtomicTransitionManager::instance = 0;

G4AtomicTransitionManager* G4AtomicTransitionManager::Instance() {
  if (instance = 0) {

    instance = new G4AtomicTransitionManager;
  }
  return instance;
}

//....................000000000000000...................000000000000............00000000000000

// this function read from EADL the binding energies of any subshell of the atom, and creates 
// the objects "shells" to be put in a map where to every Z is associated the vector of the corresponding shells

void G4AtomicTransitionManager::BuildTable() {

  G4int Z; 

  G4int dataNum = 2;
  
  G4SecondLevel* theBindingEnergyTable = util.BuildSecondLevelTables(0,dataNum,"fluor/binding");

  for (Z = 1; Z<=99; Z++) {

 
    dataNum=3;
    G4SecondLevel* theAtomShellFL = util.BuildSecondLevelTables(Z, dataNum, "fluor/fl-tr-pr-");


    //Remember first position in vector is position 0

    G4FirstLevel* theBindEnVec = (*theBindingEnergyTable)[Z-1];
    
    G4std::vector<G4AtomicShell> zShellVec;

    G4std::vector<G4int> finalTransitionSubShellIdVect(0);
    G4std::vector<G4double> transitionProbabilities(0);
    G4std::vector<G4double> transitionEnergies(0);
    
    
    for(G4int subShellId = 0; subShellId < (*theBindEnVec)[0]->size(); subShellId++) {

      G4int shellId = (G4int) (*(*theBindEnVec)[0])[subShellId];

      G4double theBindEnergy=((*(*theBindEnVec)[1])[subShellId])*MeV;

      G4FirstLevel* theSubShellFL = (*theAtomShellFL)[subShellId];


      // Fluorescence data r present only for Z > 5

      G4AtomicShell subShell;

      if(Z > 5){
      

	// theSubShellFL is a table that has, 4 each shell (in wich ther is a vacancy) of the atom three columns: 
	// the subshell availables, the probability that this transition occurs, 
	// and the energy of the gamma emitted

	const G4int SubShellCol = 0, ProbCol = 1, EnergyCol = 2;

	// Transnum is the index of the loop and of the table of transition. it starts from 1
	// because the first element of the data table is the primary shell id number and not a
	// transition probability.

	for(G4int TransNum = 1; TransNum < (*theSubShellFL)[ProbCol]->size(); TransNum++) {

	  G4int finalTrSubShell = (*(*theSubShellFL)[SubShellCol])[TransNum];
	  G4double transProb = (G4double) (*(*theSubShellFL)[ProbCol])[TransNum];
	  G4double gammaTransEnergy = (G4double) ((*(*theSubShellFL)[EnergyCol])[TransNum])*MeV;
	
	  // finTranSubShellId is a vector where there r the possible destination subshells' Ids
	  // transProbabilities is a vector where is stored the probabilities of the above destinations occurs
	  // transEnergies is a vector containing the energy of the gamma emitted by above transitions.

	  finalTransitionSubShellIdVect.push_back(finalTrSubShell); 
	  transitionProbabilities.push_back(transProb); 
	  transitionEnergies.push_back(gammaTransEnergy);
	}


	// Now I have all the data needed to create a subshell object, and i do it:

	// G4AtomicShell subShell = new G4AtomicShell(shellId,theBindEnergy,finalTransitionSubShellIdVect,
	//					   transitionProbabilities,transitionEnergies);
      }
       
    
      // describing what happens if Z<=5: subshell is created but data regarding transitions are null
      /*
      else {

	finalTransitionSubShellIdVect = 0; 
	transitionProbabilities = 0.0; 
	transitionEnergies = 0.0;

	// G4AtomicShell subShell = new G4AtomicShell(shellId,theBindEnergy,finalTransitionSubShellIdVect,
	// transitionProbabilities,transitionEnergies);
	}*/

      subShell =  G4AtomicShell(shellId,theBindEnergy,finalTransitionSubShellIdVect,
				transitionProbabilities,transitionEnergies);

      // a shell vector for the atom Z has been created: now I put it to the map of the atom's shells

      zShellVec.push_back(subShell);
    }
    /*
    for (G4int sShell = 0; sShell< zShellVec.size(); sShell++) {

      shellTable[Z][sShell] = zShellVec[sShell];
    }
    */


	     
    shellTable[Z].operator=(zShellVec); 
  }

}




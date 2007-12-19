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
#include "AIDA/AIDA.h"

int main(int argc, char* argv[]){
  
  
  G4int Z = 0;
  G4int shellId = -1;
  G4int shellIndex = -1;
  G4double totFluoProb;


  AIDA::ITree* tree;
  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITupleFactory* tupleFactory;
  AIDA::ITuple* tupleFluoProb;

  analysisFactory = AIDA_createAnalysisFactory();
  AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
  tree = treeFactory->create("Probabilities.xml","xml",false,true);
  tupleFactory = analysisFactory->createTupleFactory(*tree);
  // Book tuple column names
  std::vector<std::string> columnNames;
  // Book tuple column types
  std::vector<std::string> columnTypes;

      columnNames.push_back("AtomicNumber");
      columnNames.push_back("ShellId");
      columnNames.push_back("Probability");
      
      columnTypes.push_back("int");
      columnTypes.push_back("int");
      columnTypes.push_back("double");
      tupleFluoProb = tupleFactory->create("10", "Totale shell probabilities", columnNames, columnTypes, "");
      assert(tupleFluoProb);
  
  //  G4cout <<argc<<G4endl;
  //  G4cout <<argv[1]<<G4endl;
  //  G4cout <<argv[2]<<G4endl;

  if (argv[1]) {Z= atoi(argv[1]);}
  else {
    G4cout << "Enter Z " << G4endl;
    G4cin >> Z;
  }
  if ( argc==3 ) {shellId = atoi(argv[2]);}
  else if (argc != 2){
    G4cout <<" Select the Id of the Fluo vacancy "<<G4endl;
    G4cin>> shellId; 
  }
  
  G4AtomicTransitionManager* transManager = G4AtomicTransitionManager::Instance();
  
  G4int zStart, zEnd;
  G4bool totalFlag;
  if (Z == 0) {
    zStart = 6;
    zEnd = 100;
    totalFlag = true;
  }
  else { 
    zStart = Z;
    zEnd = Z+1;
  }
  for (Z=zStart; Z<zEnd; Z++){
    
    
    G4int shellNumber = transManager->NumberOfReachableShells(Z);



    for (shellIndex=0; shellIndex<shellNumber; shellIndex++) {

      if (shellId == -1) {
	shellIndex = -1;
	break;
      }
      if (transManager->Shell(Z,shellIndex)->ShellId() == shellId)  {
	G4cout << "ShellId "<< shellId <<" matches the index: " << shellIndex << G4endl;
	break;
      }
      
      else if (shellIndex == shellNumber-1 ) {
	G4cout << "ShellId " << shellId << " not available. " << G4endl;
	return 0;
      }
     }
    
    
    G4int shellStart, shellEnd;
    if (shellIndex == -1 ) {
      shellStart = 0;
      shellEnd = shellNumber;
    }
    else { 
      shellStart = shellIndex;
      shellEnd = shellIndex+1;
    }
    
    
    
    for (shellIndex=shellStart; shellIndex<shellEnd; shellIndex++) { 
      
      // G4cout << "Testing G4FluoTransition "<<G4endl;
      //  std::vector<G4double> transEnergies = transManager->ReachableShell(Z,shellIndex)->TransitionEnergies();
      //  std::vector<G4int> transIds = transManager->ReachableShell(Z,shellIndex)->OriginatingShellIds();
      //std::vector<G4double> transProbs = transManager->ReachableShell(Z,shellIndex)->TransitionProbabilities();
      G4int currentShellId = transManager->Shell(Z,shellIndex)->ShellId();
      if (shellId == -1) {
	G4cout << "ShellId: "<< currentShellId;
      }
      totFluoProb = transManager->TotalRadiativeTransitionProbability(Z,shellIndex);

      G4cout << "Total Probability that a radiative transition occurs: " <<
	totFluoProb <<G4endl;
      if (totalFlag) {
	tupleFluoProb->fill(0,Z);
	tupleFluoProb->fill(1,currentShellId);
	tupleFluoProb->fill(2,totFluoProb);
	tupleFluoProb->addRow();
      }
    }
  }
  tree->commit(); // Write histos in file. 
  tree->close();
}

  //  G4cout << "Number of shells: " << transManager->NumberOfShells(Z)<<G4endl;

  //  G4cout<< "Number of Fluo vacancies: "<<transManager->NumberOfReachableShells(Z)<<G4endl;

  //  G4cout << "Testing G4AtomicShell" << G4endl;




  //  G4cout << "Select the index of the subshell whose binding energy you need: " << G4endl;
  //  G4cin >> subShell;



//  G4cout << "SubShell binding energy: " << transManager->Shell(Z,subShell)->BindingEnergy() << G4endl;
  
//  G4cout << "Testing G4AtomicShell" << G4endl;
//  G4int shellIndex;



//  for (G4int trans=0; trans<transIds.size(); trans++) {
//
//    G4cout << "The transition starts from the shell: " << transIds[trans] << G4endl;
//    G4cout << "The transition starts from the shell: " << 
//
//
//      transManager->ReachableShell(Z,shellIndex)->OriginatingShellId(trans) << G4endl;
//
//
//    G4cout << " Transition energy: " << transEnergies[trans] << G4endl;
//    G4cout << "Transition energy: " << 
//      transManager->ReachableShell(Z,shellIndex)->TransitionEnergy(trans) << G4endl;
//
//    G4cout << "Transition probability: " << transProbs[trans] << G4endl;
//    G4cout << "Transition probability: " << 
//      transManager->ReachableShell(Z,shellIndex)->TransitionProbability(trans) << G4endl;
//  }

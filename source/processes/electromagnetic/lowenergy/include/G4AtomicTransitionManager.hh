// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4TransitionManager.hh
//
//      Author:        Alfonso Mantero (alf@mailandnews.com)
// 
//      Creation date: 1 May 2001
//
// -------------------------------------------------------------------

#ifndef G4AtomicTransition_h
#define G4AtomicTransition_h 1

#include "G4LowEnergyUtilities.hh"
#include "G4AtomicShell.hh"
#include "globals.hh"
#include "g4std/algorithm"
#include "g4std/map"
#include "g4std/vector"

//this class is a SINGLETON
class G4AtomicTransitionManager {

public: 

  static G4AtomicTransitionManager* Instance();
  //the only way to get an instance of this class is to call the 
  //function Instance() 
  const G4AtomicShell* Shell(G4int z, G4int shellIdentifier);
  //z is the atomic number of the element, shellIdentifier is the 
  //name (in EADL) of the starting shell for the transition
  G4int NumberOfShells(G4int z);
  //this function returns the number of shells of the element
  //whose atomic number is z
  G4double TotalRadiativeTransitionProbability(G4int z, G4int shellId);
  //gives the sum of the probabilities of radiative transition from the
  //shell named shellId of an atom og the element z
  G4double TotalNonRadiativeTransitionProbability(G4int z, G4int shellId);
 //gives the sum of the probabilities of non radiative transition from the
  //shell named shellId of an atom og the element z

protected:

  G4AtomicTransitionManager();
  ~G4AtomicTransitionManager();

private:
  G4LowEnergyUtilities util;
  static G4AtomicTransitionManager* instance;
  
  G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> > shellTable;
  // the first element of the map is the atomic number z.
  //the second element is a vector of shells.
  
  G4std::vector< G4int >ZVector;

  //this vector contains all the atomic numbers of the elements contained
  //in the MaterialTable. Their order is the same in wich they appeare in the
  //table

  G4int tableLimit;
  //EADL contains the data for radiative transition only for atoms with Z>5.
  //tableLimit is initialized to 5 in the constructor of this class
 
   void BuildZVec();

  G4SecondLevel*BuildBindingEnergiesTable();
  G4ThirdLevel*BuildTransitionTable();
  //these two functions read the files of EADL's data and build some table
};

#endif

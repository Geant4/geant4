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
// $Id: G4AtomicTransitionManager.hh,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Authors: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//          Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// ?????      Created
//
// -------------------------------------------------------------------

#ifndef G4AtomicTransitionManager_h
#define G4AtomicTransitionManager_h 1

#include "G4ShellData.hh"
#include "G4FluoData.hh"
#include "G4AtomicTransition.hh"
#include "G4AtomicShell.hh"
#include "globals.hh"
#include "g4std/algorithm"
#include "g4std/map"
#include "g4std/vector"

//this class is a SINGLETON
class G4AtomicTransitionManager {

public: 

  static G4AtomicTransitionManager* Instance();
  // The only way to get an instance of this class is to call the 
  // function Instance() 

  const G4AtomicShell* Shell(G4int Z, size_t shellIndex);
  //Z is the atomic number of the element, shellIdentifier is the 
  //index (in EADL) of the shell 

  const G4AtomicTransition* ReachableShell(G4int Z, size_t shellIndex);
  //Z is the atomic number of the element, shellIdentifier is the 
  //index (in EADL) of the final shell for the transition 

  G4int NumberOfShells(G4int Z);
  //this function returns the number of shells of the element
  //whose atomic number is Z

  G4int NumberOfReachableShells(G4int Z);

  G4double TotalRadiativeTransitionProbability(G4int Z, size_t shellIndex);
  //gives the sum of the probabilities of radiative transition towards the
  //shell whose index is shellId
 
  G4double TotalNonRadiativeTransitionProbability(G4int Z, size_t shellIndex);
  //gives the sum of the probabilities of non radiative transition from the
  //shell whose index is shellId 


protected:

  G4AtomicTransitionManager(G4int minZ = 1, G4int maxZ = 99, 
			    G4int limitInfTable = 6, G4int limitSupTable=100 );
  ~G4AtomicTransitionManager();

private:
 
  static G4AtomicTransitionManager* instance;
  
  G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> > shellTable;
  // the first element of the map is the atomic number z.
  //the second element is a vector of shells.

  G4std::map<G4int,G4std::vector<G4AtomicTransition*>,std::less<G4int> > transitionTable;
  
  G4int zMin;
  G4int zMax;
  G4int infTableLimit;
  G4int supTableLimit;
  //EADL contains the data for radiative transition only for atoms with Z>5.
  //tableLimit is initialized to 5 in the constructor of this class
 
};

#endif

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
//
//
// Authors: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//          Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  
//  16 Sept 2001 EG  Modified according to a design iteration in the 
//                   LowEnergy category
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics: create or fills and manages G4RDAtomicShell, 
// G4RDFluoTransition, G4RDAugerTransition objects.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDAtomicTransitionManager_h
#define G4RDAtomicTransitionManager_h 1

#include "G4RDShellData.hh"
#include "G4RDFluoData.hh"
#include "G4RDAugerData.hh"
#include "G4RDFluoTransition.hh"
#include "G4RDAtomicShell.hh"
// #include "g4std/map"
#include <vector>
#include "globals.hh"

// This class is a singleton
class G4RDAtomicTransitionManager {

public: 

  // The only way to get an instance of this class is to call the 
  // function Instance() 
  static G4RDAtomicTransitionManager* Instance();
 
  // Z is the atomic number of the element, shellIndex is the 
  // index (in EADL) of the shell
  G4RDAtomicShell* Shell(G4int Z, size_t shellIndex) const;
   
  // Z is the atomic number of the element, shellIndex is the 
  // index (in EADL) of the final shell for the transition
  // This function gives, upon Z and the Index of the initial shell where te vacancy is, 
  // the radiative transition that can happen (originating shell, energy, probability)
  const G4RDFluoTransition* ReachableShell(G4int Z, size_t shellIndex) const ;

  // This function gives, upon Z and the Index of the initial shell where te vacancy is, 
  // the NON-radiative transition that can happen with originating shell for the transition, and the
  // data for the possible auger electrons emitted (originating vacancy, energy amnd probability)
  
  const G4RDAugerTransition* ReachableAugerShell(G4int Z, G4int shellIndex) const ;
  
  // This function returns the number of shells of the element
  // whose atomic number is Z
  G4int NumberOfShells(G4int Z) const;
  
  // This function returns the number of those shells of the element
  // whose atomic number is Z which are reachable through a radiative
  // transition
  
  // This function returns the number of possible radiative transitions for the atom with atomic number Z
  // i.e. the number of shell in wich a vacancy can be filled with a radiative transition 
  
  G4int NumberOfReachableShells(G4int Z)const ;

// This function returns the number of possible NON-radiative transitions for the atom with atomic number Z
// i.e. the number of shell in wich a vacancy can be filled by a NON-radiative transition 

  G4int NumberOfReachableAugerShells(G4int Z)const ;

  // Gives the sum of the probabilities of radiative transition towards the
  // shell whose index is shellIndex
  G4double TotalRadiativeTransitionProbability(G4int Z, size_t shellIndex);
  
  // Gives the sum of the probabilities of non radiative transition from the
  // shell whose index is shellIndex
  G4double TotalNonRadiativeTransitionProbability(G4int Z, size_t shellIndex);
   
protected:

  G4RDAtomicTransitionManager(G4int minZ = 1, G4int maxZ = 100, 
			    G4int limitInfTable = 6, G4int limitSupTable=100 );
  ~G4RDAtomicTransitionManager();

private:
  // Hide copy constructor and assignment operator 
  G4RDAtomicTransitionManager& operator=(const G4RDAtomicTransitionManager& right);
  G4RDAtomicTransitionManager(const G4RDAtomicTransitionManager&);
 
  static G4RDAtomicTransitionManager* instance;
  
  // the first element of the map is the atomic number Z.
  // the second element is a vector of G4RDAtomicShell*.
  std::map<G4int,std::vector<G4RDAtomicShell*>,std::less<G4int> > shellTable;
  
  // the first element of the map is the atomic number Z.
  // the second element is a vector of G4AtomicTransition*.
  std::map<G4int,std::vector<G4RDFluoTransition*>,std::less<G4int> > transitionTable;
 
  // since Augereffect data r stored as a table in G4RDAugerData, we have here a pointer to an element of that class itself.

  G4RDAugerData* augerData;

  // Minimum and maximum Z in EADL table containing identities and binding
  // energies of shells
  G4int zMin;
  G4int zMax;
  
  // Minimum and maximum Z in EADL table containing identities, transition 
  // energies and transition probabilities of shells
  G4int infTableLimit;
  G4int supTableLimit;
 
 
};

#endif

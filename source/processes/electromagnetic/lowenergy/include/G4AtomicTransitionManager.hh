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
// Low Energy Electromagnetic Physics: create or fills and manages G4AtomicShell, 
// G4FluoTransition, G4AugerTransition objects.
// -------------------------------------------------------------------

#ifndef G4AtomicTransitionManager_h
#define G4AtomicTransitionManager_h 1

#include "G4ShellData.hh"
#include "G4FluoTransition.hh"
#include "G4AugerTransition.hh"
#include "G4AtomicShell.hh"
#include <vector>
#include "globals.hh"

class G4AugerData;

// This class is a singleton
class G4AtomicTransitionManager {

public: 
  /// The only way to get an instance of this class is to call the 
  /// function Instance() 
  static G4AtomicTransitionManager* Instance();

  /// needs to be called once from other code before start of run
  void Initialise();
 
  /// Z is the atomic number of the element, shellIndex is the 
  /// index (in EADL) of the shell
  G4AtomicShell* Shell(G4int Z, size_t shellIndex) const;
   
  /// Z is the atomic number of the element, shellIndex is the 
  /// index (in EADL) of the final shell for the transition
  /// This function gives, upon Z and the Index of the initial shell where 
  /// the vacancy is, the radiative transition that can happen (originating 
  /// shell, energy, probability)
  const G4FluoTransition* ReachableShell(G4int Z, size_t shellIndex) const;

  /// This function gives, upon Z and the Index of the initial shell where 
  /// the vacancy is, the NON-radiative transition that can happen with 
  /// originating shell for the transition, and the data for the possible 
  /// auger electrons emitted (originating vacancy, energy amnd probability)  
  const G4AugerTransition* ReachableAugerShell(G4int Z, G4int shellIndex) const;
  
  /// This function returns the number of shells of the element
  /// whose atomic number is Z
  G4int NumberOfShells(G4int Z) const;
  
  /// This function returns the number of those shells of the element
  /// whose atomic number is Z which are reachable through a radiative
  /// transition
  G4int NumberOfReachableShells(G4int Z) const;

  /// This function returns the number of possible NON-radiative transitions 
  /// for the atom with atomic number Z i.e. the number of shell in wich 
  /// a vacancy can be filled by a NON-radiative transition 
  G4int NumberOfReachableAugerShells(G4int Z) const;

  /// Gives the sum of the probabilities of radiative transition towards the
  /// shell whose index is shellIndex
  G4double 
  TotalRadiativeTransitionProbability(G4int Z, size_t shellIndex) const;
  
  /// Gives the sum of the probabilities of non radiative transition from the
  /// shell whose index is shellIndex
  G4double 
  TotalNonRadiativeTransitionProbability(G4int Z, size_t shellIndex) const; 

  /// Verbosity control
  void SetVerboseLevel(G4int vl) {verboseLevel = vl;};
  G4int GetVerboseLevel(){return verboseLevel;};

private:
  explicit G4AtomicTransitionManager();

  ~G4AtomicTransitionManager();

  // Hide copy constructor and assignment operator 
  G4AtomicTransitionManager& operator=(const G4AtomicTransitionManager& right);
  G4AtomicTransitionManager(const G4AtomicTransitionManager&);
 
  static G4AtomicTransitionManager* instance;
  // since Augereffect data r stored as a table in G4AugerData, we have 
  // here a pointer to an element of that class itself.
  G4AugerData* augerData;

  // the first element of the map is the atomic number Z.
  // the second element is a vector of G4AtomicShell*.
  std::map<G4int,std::vector<G4AtomicShell*>,std::less<G4int> > shellTable;
  
  // the first element of the map is the atomic number Z.
  // the second element is a vector of G4AtomicTransition*.
  std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> > transitionTable;
 
  // Minimum and maximum Z in EADL table containing identities and binding
  // energies of shells
  G4int zMin = 1; 
  G4int zMax = 104;
  
  // Minimum and maximum Z in EADL table containing identities, transition 
  // energies and transition probabilities of shells
  G4int infTableLimit = 6;
  G4int supTableLimit = 104;
  G4int verboseLevel;
  G4bool isInitialized;
};

#endif

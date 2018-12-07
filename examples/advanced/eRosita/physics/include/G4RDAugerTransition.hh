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
// 
// Author: Alfonso Mantero (Alfosno.Mantero@ge.infn.it)
//
// History:
// -----------
// 1 Jun 2002: first Commited to CVS
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// This Class stores all the information of auger effect relative 
// to one main vacancy in a atom, like possible auger emission with 
// relative probabilites, originating shell's Ids, probabilities of 
// transition and auger electron energies. 
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDAugerTransition_h 
#define G4RDAugerTransition_h 1

#include "G4DataVector.hh"
#include "globals.hh"
#include <vector>
#include <map>

class G4RDAugerTransition {

public:

  G4RDAugerTransition(G4int finalShell, std::vector<G4int> transIds,
		    const std::map<G4int, std::vector<G4int>, std::less<G4int> >* idMap,
		    const std::map<G4int, G4DataVector, std::less<G4int> >* energyMap,
		    const std::map<G4int, G4DataVector, std::less<G4int> >* probabilityMap);

  ~G4RDAugerTransition();
  
// All the data stored and provided by this class are relative to a
// given vacancy, whose identity is provided by the FinalShellId() method,
// in an atom of a given material

// Returns the ids of the shells from wich an auger electron culd came from, given the shell
// from wich the transition electron comes from.

  const std::vector<G4int>* AugerOriginatingShellIds(G4int startShellId) const;

// Returns the ids of the shells from wich an electron cuuld fill the vacancy in finalShellId

  const std::vector<G4int>* TransitionOriginatingShellIds() const;

// Returns the energiess of the possible auger electrons, given th shell
// from wich the transition electron comes from.

  const G4DataVector* AugerTransitionEnergies(G4int startShellId) const;

// Returns the emission probabilities of the auger electrons, given th shell
// from wich the transition electron comes from.

  const G4DataVector* AugerTransitionProbabilities(G4int startShellId) const;

// returns the id of the shell in wich the transition electron arrives

  G4int FinalShellId() const;

// Returns the id of the shell from wich come the auger electron , given the shell
// from wich the transition electron comes from and the index number.

  G4int AugerOriginatingShellId(G4int index, G4int startShellId) const;

// Returns the energy of the auger electron, given the shell
// from wich the transition electron comes from and the index number.

  G4double AugerTransitionEnergy(G4int index, G4int startShellId) const;

// Returns the probability of the auger emission, given the shell
// from wich the transition electron comes from and the index number.

  G4double AugerTransitionProbability(G4int index, G4int startShellId) const;

// Returns the id of the shell form wich the transition electron come from

  G4int TransitionOriginatingShellId(G4int index) const;


private:

  G4int finalShellId;
  std::map<G4int,std::vector<G4int>,std::less<G4int> >  augerOriginatingShellIdsMap;
  std::map<G4int,G4DataVector,std::less<G4int> >  augerTransitionEnergiesMap;
  std::map<G4int,G4DataVector,std::less<G4int> >  augerTransitionProbabilitiesMap;
  std::vector<G4int> transitionOriginatingShellIds;
  
};

#endif



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
// $Id: G4AugerTransition.hh v0.1
// 
// 
// Author: Alfonso Mantero (Alfosno.Mantero@ge.infn.it)
//
// History:
// -----------
// 6 Mar 2002: first implementation
//
// -------------------------------------------------------------------


#ifndef G4AugerTransition_h 
#define G4AugerTransition_h 1

#include "G4DataVector.hh"
#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"

class G4AugerTransition {

public:

  G4AugerTransition(G4int finalShell, G4std::vector<G4int> transIds,
		    const G4std::map<G4int, G4std::vector<G4int>, G4std::less<G4int> >* idMap,
		    const G4std::map<G4int, G4DataVector, G4std::less<G4int> >* energyMap,
		    const G4std::map<G4int, G4DataVector, G4std::less<G4int> >* probabilityMap);

  ~G4AugerTransition();
  
// All the data stored and provided by this class are relative to a
// given vacancy, whose identity is provided by the FinalShellId() method,
// in an atom of a given material

// Returns the ids of the shells from wich an auger electron culd came from, given the shell
// from wich the transition electron comes from.

  const G4std::vector<G4int>* AugerOriginatingShellIds(G4int startShellId) const;

// Returns the ids of the shells from wich an electron cuuld fill the vacancy in finalShellId

  const G4std::vector<G4int>* TransitionOriginatingShellIds() const;

// Returns the energiess of the possible auger electrons, given th shell
// from wich the transition electron comes from.

  const G4DataVector* AugerTransitionEnergies(G4int startShellId) const;

// Returns the emission probabilities of the auger electrons, given th shell
// from wich the transition electron comes from.

  const G4DataVector* AugerTransitionProbabilities(G4int startShellId) const;

// returns the id of the shell in wich the transition electron arrives

  const G4int FinalShellId() const;

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
  G4std::map<G4int,G4std::vector<G4int>,G4std::less<G4int> >  augerOriginatingShellIdsMap;
  G4std::map<G4int,G4DataVector,G4std::less<G4int> >  augerTransitionEnergiesMap;
  G4std::map<G4int,G4DataVector,G4std::less<G4int> >  augerTransitionProbabilitiesMap;
  G4std::vector<G4int> transitionOriginatingShellIds;
  
};

#endif



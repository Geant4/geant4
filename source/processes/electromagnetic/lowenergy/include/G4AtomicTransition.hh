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
// $Id: G4AtomicTransition.hh,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  16 Sept 2001 Modofied according to a design iteration in the 
//              LowEnergy category
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, a data container
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------


#ifndef G4AtomicTransition_h 
#define G4AtomicTransition_h 1

#include "G4DataVector.hh"
#include "globals.hh"
#include "g4std/vector"

class G4AtomicTransition {

public:

  G4AtomicTransition(G4int finalShell,
		     const G4std::vector<G4int>& ids,
		     const G4DataVector& energies,
		     const G4DataVector& probabilities);
  ~G4AtomicTransition();
  
  // All the data stored and provided by this class are relative to a
  // given vacancy, whose identity is provided by the FinalShellId() method,
  // in an atom of a given material

  // Returns the identities of the originating shells for the transitions 
  const G4std::vector<G4int>& OriginatingShellIds() const;
  
  // Return the energies of the transitions
  const G4DataVector& TransitionEnergies() const;

  // Return the probabilities of the transitions
  const G4DataVector& TransitionProbabilities() const;
  
  // Return the identity if the vacancy
  const G4int FinalShellId() const;

  // Given the index of the originating shells returns its identity
  G4int OriginatingShellId(G4int index) const;

  // Given the index of the originating shells returns the energy
  // of the transition starting from it
  G4double TransitionEnergy(G4int index) const;

  // Given the index of the originating shells returns the probability
  // of the transition starting from it
  G4double TransitionProbability(G4int index) const;

private:

  G4int finalShellId;
  G4std::vector<G4int> originatingShellIds;
  G4DataVector transitionEnergies;
  G4DataVector transitionProbabilities;
  
};

#endif


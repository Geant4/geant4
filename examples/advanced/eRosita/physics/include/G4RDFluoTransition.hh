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


#ifndef G4RDFluoTransition_h 
#define G4RDFluoTransition_h 1

#include "G4DataVector.hh"
#include "globals.hh"
#include <vector>

class G4RDFluoTransition {

public:

  G4RDFluoTransition(G4int,const std::vector<G4int>&,const G4DataVector&,
		     const G4DataVector&);

  ~G4RDFluoTransition();
  
  // All the data stored and provided by this class are relative to a
  // given vacancy, whose identity is provided by the FinalShellId() method,
  // in an atom of a given material

  // Returns the identities of the originating shells for the transitions 
  const std::vector<G4int>& OriginatingShellIds() const;
  
  // Return the energies of the transitions
  const G4DataVector& TransitionEnergies() const;

  // Return the probabilities of the transitions
  const G4DataVector& TransitionProbabilities() const;
  
  // Return the identity if the vacancy
  G4int FinalShellId() const;

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
  std::vector<G4int> originatingShellIds;
  G4DataVector transitionEnergies;
  G4DataVector transitionProbabilities;
  
};

#endif


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
// $Id$
//
// Multibody "phase space" generator, which provides interface to
// algorithms for sampling.  Momentum vectors are generated in the
// center-of-mass frame of the decay, and returned in a user-supplied
// buffer.  
//
// Sampling algorithm is specified via constructor argument, either by
// code for centrally supplied algorithms (Kopylov, GENBOD, or NBody),
// or by pointer (which is owned by the generator).
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#ifndef G4HadDecayGenerator_HH
#define G4HadDecayGenerator_HH 1

#include "globals.hh"
#include "G4LorentzVector.hh"
#include <vector>

class G4ParticleDefinition;
class G4VHadDecayAlgorithm;


class G4HadDecayGenerator {
public:
  // Flags to select algorithm by code in constructor
  enum Algorithm { NONE=0, Kopylov=1, GENBOD=2, NBody=3 };

  // Specify "standard" algorithm by code or by object (takes ownership)
  G4HadDecayGenerator(Algorithm alg=Kopylov, G4int verbose=0);
  G4HadDecayGenerator(G4VHadDecayAlgorithm* alg, G4int verbose=0);
  virtual ~G4HadDecayGenerator();

  // Enable (or disable if 0) diagnostic messages
  void SetVerboseLevel(G4int verbose);
  const G4String& GetAlgorithmName() const;

  // Initial state (rest mass) and list of final masses
  G4bool Generate(G4double initialMass,
		  const std::vector<G4double>& masses,
		  std::vector<G4LorentzVector>& finalState);

  // Initial state particle and list of final masses
  G4bool Generate(const G4ParticleDefinition* initialPD,
		  const std::vector<G4double>& masses,
		  std::vector<G4LorentzVector>& finalState);

  // Initial state (frame) and list of final masses
  // Final state particles will be boosted to initial-state frame
  G4bool Generate(const G4LorentzVector& initialState,
		  const std::vector<G4double>& masses,
		  std::vector<G4LorentzVector>& finalState);

protected:
  // Special case for one-body final state
  G4bool GenerateOneBody(G4double initialMass,
			 const std::vector<G4double>& masses,
			 std::vector<G4LorentzVector>& finalState) const;

  // Special function used by constructor for unrecognized algorithm code
  void ReportInvalidAlgorithm(Algorithm alg) const;
  void ReportMissingAlgorithm() const;

protected:
  // SPECIAL FUNCTION FOR SUBCLASSES: A subclass may implement a
  // collection of algorithms, to be switched on an event-by-event
  // basis.  This function allows the subclass to switch the "active"
  // algorithm before Generate() is called.
  //
  // If this function is used by the subclass, then the subclass has
  // ownership of _all_ instantiated algorithms, and should delete
  // them in its own dtor.  The subclass dtor must also call
  // UseAlgorithm(0) to set the base algorithm to a null pointer, to
  // prevent a double-delete error.
  void UseAlgorithm(G4VHadDecayAlgorithm* alg) { theAlgorithm = alg; }

  G4int verboseLevel;
  G4VHadDecayAlgorithm* theAlgorithm;
};

#endif	/* G4HadDecayGenerator_HH */

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
// Abstract base class for multibody "phase space" generators.  Subclasses
// implement a specific algorithm, such as Kopylov, GENBOD, or Makoto's
// NBody.  Subclasses are used by G4HadPhaseSpaceGenerator.
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#ifndef G4VHadDecayAlgorithm_HH
#define G4VHadDecayAlgorithm_HH 1

#include "globals.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <iosfwd>

class G4VHadDecayAlgorithm {
public:
  G4VHadDecayAlgorithm(const G4String& algName, G4int verbose=0)
    : name(algName), verboseLevel(verbose) {;}
  virtual ~G4VHadDecayAlgorithm() {;}

  // Initial state (rest mass) and list of final masses
  void Generate(G4double initialMass,
		const std::vector<G4double>& masses,
		std::vector<G4LorentzVector>& finalState);

  // Enable (or disable if 0) diagnostic messages (subclass may overload)
  virtual void  SetVerboseLevel(G4int verbose) { verboseLevel = verbose; }
  G4int GetVerboseLevel() const { return verboseLevel; }
  const G4String& GetName() const { return name; }
  
protected:
  // Subclasses MUST implement these functions
  virtual void GenerateTwoBody(G4double initialMass,
			       const std::vector<G4double>& masses,
			       std::vector<G4LorentzVector>& finalState) = 0;

  virtual void GenerateMultiBody(G4double initialMass,
				 const std::vector<G4double>& masses,
				 std::vector<G4LorentzVector>& finalState) = 0;

  // Validate kinematics (e.g., limit number of final state particles)
  // Subclasses may override or call back to this function
  virtual G4bool IsDecayAllowed(G4double initialMass,
				const std::vector<G4double>& masses) const;

  // Two-body momentum function (c.f. PDK from CERNLIB W505)
  G4double TwoBodyMomentum(G4double M0, G4double M1, G4double M2) const;

  // Convenience functions for uniform angular distributions
  G4double UniformTheta() const;
  G4double UniformPhi() const;

  // Utility to dump vector contents to line of output
  void PrintVector(const std::vector<G4double>& v, const G4String& name,
		   std::ostream& os) const;

private:
  G4String name;
  G4int verboseLevel;
};

#endif	/* G4VHadDecayAlgorithm_HH */

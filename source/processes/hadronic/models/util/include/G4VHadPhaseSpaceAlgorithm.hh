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
// Abstract base class for multibody uniform phase space generators.
// Subclasses implement a specific algorithm, such as Kopylov, GENBOD,
// or Makoto's NBody.  Subclasses are used by G4HadDecayGenerator.
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#ifndef G4VHadPhaseSpaceAlgorithm_HH
#define G4VHadPhaseSpaceAlgorithm_HH 1

#include "globals.hh"
#include "G4VHadDecayAlgorithm.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <iosfwd>

class G4VHadPhaseSpaceAlgorithm : public G4VHadDecayAlgorithm {
public:
  G4VHadPhaseSpaceAlgorithm(const G4String& algName, G4int verbose=0)
    : G4VHadDecayAlgorithm(algName, verbose) {;}
  virtual ~G4VHadPhaseSpaceAlgorithm() {;}

protected:
  // Multi-body function remains pure-virtual from base class

  // Two-body uniform distribution is common to all algorithms
  virtual void GenerateTwoBody(G4double initialMass,
			       const std::vector<G4double>& masses,
			       std::vector<G4LorentzVector>& finalState);

  // Sample spherical distribution
  G4ThreeVector UniformVector(G4double mag=1.) const;
};

#endif	/* G4VHadPhaseSpaceAlgorithm_HH */

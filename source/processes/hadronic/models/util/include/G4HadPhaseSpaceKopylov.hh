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
// Multibody "phase space" generator using Kopylov's algorithm
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#ifndef G4HadPhaseSpaceKopylov_HH
#define G4HadPhaseSpaceKopylov_HH 1

#include "G4VHadPhaseSpaceAlgorithm.hh"


class G4HadPhaseSpaceKopylov : public G4VHadPhaseSpaceAlgorithm {
public:
  G4HadPhaseSpaceKopylov(G4int verbose=0)
    : G4VHadPhaseSpaceAlgorithm("G4HadPhaseSpaceKopylov",verbose) {;}
  virtual ~G4HadPhaseSpaceKopylov() {;}

protected:
  virtual void GenerateMultiBody(G4double initialMass,
				 const std::vector<G4double>& masses,
				 std::vector<G4LorentzVector>& finalState);

protected:
  G4double BetaKopylov(G4int K) const;
};

#endif	/* G4HadPhaseSpaceKopylov_HH */

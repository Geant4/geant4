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
// Hadronic Process: Phase space decay for the Fermi BreakUp model
// by V. Lara
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko: 
//          - IsotropicVector is inlined
//          - Momentum computation return zero or positive value
//          - DumpProblem method is added providing more information
//          - Reduced usage of exotic std functions  
//

#ifndef G4FermiPhaseSpaceDecay_hh
#define G4FermiPhaseSpaceDecay_hh 1

#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <vector>

class G4Pow;

class G4FermiPhaseSpaceDecay
{
public:

  G4FermiPhaseSpaceDecay();
  ~G4FermiPhaseSpaceDecay();
  
  std::vector<G4LorentzVector*>* Decay(G4double parent_mass, 
	const std::vector<G4double>& fragment_masses) const;

private:

  inline G4double PtwoBody(G4double E, G4double P1, G4double P2) const;
  
  G4double BetaKopylov(G4int, CLHEP::HepRandomEngine*) const; 

  std::vector<G4LorentzVector*> * 
  KopylovNBodyDecay(G4double, const std::vector<G4double>&) const;

  G4FermiPhaseSpaceDecay(const G4FermiPhaseSpaceDecay&) = delete;
  const G4FermiPhaseSpaceDecay & operator=
  (const G4FermiPhaseSpaceDecay &) = delete; 
  G4bool operator==(const G4FermiPhaseSpaceDecay&) = delete;
  G4bool operator!=(const G4FermiPhaseSpaceDecay&) = delete;

  G4Pow* g4calc;
};

inline G4double 
G4FermiPhaseSpaceDecay::PtwoBody(G4double E, G4double P1, G4double P2) const
{
  G4double P = (E+P1+P2)*(E+P1-P2)*(E-P1+P2)*(E-P1-P2)/(4.0*E*E);
  return (P>0.0) ? std::sqrt(P) : 0.0; 
}

#endif

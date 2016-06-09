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
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4FermiPhaseSpaceDecay_hh
#define G4FermiPhaseSpaceDecay_hh


#include "G4LorentzVector.hh"
#include "G4ParticleMomentum.hh"
#include <CLHEP/Random/RandGamma.h>

#include <vector>
#include <cmath>

class G4FermiPhaseSpaceDecay
{
public:
  inline G4FermiPhaseSpaceDecay();
  inline ~G4FermiPhaseSpaceDecay();
  
  inline std::vector<G4LorentzVector*> * 
  Decay(const G4double,  const std::vector<G4double>&) const;

private:
  inline G4FermiPhaseSpaceDecay(const G4FermiPhaseSpaceDecay&);
  inline const G4FermiPhaseSpaceDecay & operator=(const G4FermiPhaseSpaceDecay &); 
  inline G4bool operator==(const G4FermiPhaseSpaceDecay&);
  inline G4bool operator!=(const G4FermiPhaseSpaceDecay&);

  inline G4double PtwoBody(G4double E, G4double P1, G4double P2) const;
  
  G4ParticleMomentum IsotropicVector(const G4double Magnitude = 1.0) const;

  inline G4double BetaKopylov(const G4int) const; 

  std::vector<G4LorentzVector*> * 
  TwoBodyDecay(const G4double, const std::vector<G4double>&) const;

  std::vector<G4LorentzVector*> * 
  NBodyDecay(const G4double, const std::vector<G4double>&) const;

  std::vector<G4LorentzVector*> * 
  KopylovNBodyDecay(const G4double, const std::vector<G4double>&) const;
};

#include "G4FermiPhaseSpaceDecay.icc"

#endif

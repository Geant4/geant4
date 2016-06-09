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
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4FermiPhaseSpaceDecay_hh
#define G4FermiPhaseSpaceDecay_hh


#include "G4LorentzVector.hh"
#include "G4ParticleMomentum.hh"
#include <CLHEP/Random/RandGamma.h>

#include <vector>


#define G4RandGamma  RandGamma

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

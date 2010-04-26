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
#include "Randomize.hh"
#include "G4Pow.hh"
#include <CLHEP/Random/RandGamma.h>

#include <vector>
#include <cmath>

class G4FermiPhaseSpaceDecay
{
public:

  G4FermiPhaseSpaceDecay();
  ~G4FermiPhaseSpaceDecay();
  
  inline std::vector<G4LorentzVector*> * 
  Decay(const G4double,  const std::vector<G4double>&) const;

private:

  G4FermiPhaseSpaceDecay(const G4FermiPhaseSpaceDecay&);
  const G4FermiPhaseSpaceDecay & operator=(const G4FermiPhaseSpaceDecay &); 
  G4bool operator==(const G4FermiPhaseSpaceDecay&);
  G4bool operator!=(const G4FermiPhaseSpaceDecay&);

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

inline G4double 
G4FermiPhaseSpaceDecay::PtwoBody(G4double E, G4double P1, G4double P2) const
{
  G4double res = -1.0;
  G4double P = (E+P1+P2)*(E+P1-P2)*(E-P1+P2)*(E-P1-P2)/(4.0*E*E);
  if (P>0.0) { res = std::sqrt(P); }
  return res;
}

inline std::vector<G4LorentzVector*> * G4FermiPhaseSpaceDecay::
Decay(const G4double parent_mass, const std::vector<G4double>& fragment_masses) const
{
  return KopylovNBodyDecay(parent_mass,fragment_masses);
}

inline G4double 
G4FermiPhaseSpaceDecay::BetaKopylov(const G4int K) const
{
  //JMQ 250410 old algorithm has been commented
  // Notice that alpha > beta always
  // const G4double beta = 1.5;
  // G4double alpha = 1.5*(K-1);
  // G4double Y1 = CLHEP::RandGamma::shoot(alpha,1);
  // G4double Y2 = CLHEP::RandGamma::shoot(beta,1);
  
  // return Y1/(Y1+Y2);

  G4Pow* g4pow = G4Pow::GetInstance();
  G4int N = 3*K - 5;
  G4double xN = G4double(N);
  G4double F;
  //G4double Fmax = std::pow((3.*K-5.)/(3.*K-4.),(3.*K-5.)/2.)*std::sqrt(1-((3.*K-5.)/(3.*K-4.))); 
  // VI variant
  G4double Fmax = std::sqrt(g4pow->powZ(N, xN/(xN + 1))/(xN + 1)); 
  G4double chi;
  do 
    {
      chi = G4UniformRand();
      F = std::sqrt(g4pow->powZ(N, chi)*(1-chi));      
    } while ( Fmax*G4UniformRand() > F);  
  return chi;

}

#endif

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
// $Id: G4ExcitationHandler.hh,v 1.13 2010-11-17 16:20:31 vnivanch Exp $
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

#include <vector>
#include <CLHEP/Units/PhysicalConstants.h>

#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4Pow.hh"

class G4FermiPhaseSpaceDecay
{
public:

  G4FermiPhaseSpaceDecay();
  ~G4FermiPhaseSpaceDecay();
  
  inline std::vector<G4LorentzVector*> * 
  Decay(const G4double,  const std::vector<G4double>&) const;

private:

  inline G4double PtwoBody(G4double E, G4double P1, G4double P2) const;
  
  inline G4ThreeVector IsotropicVector(const G4double Magnitude = 1.0) const;

  inline G4double BetaKopylov(G4int) const; 

  std::vector<G4LorentzVector*> * 
  TwoBodyDecay(G4double, const std::vector<G4double>&) const;

  std::vector<G4LorentzVector*> * 
  NBodyDecay(G4double, const std::vector<G4double>&) const;

  std::vector<G4LorentzVector*> * 
  KopylovNBodyDecay(G4double, const std::vector<G4double>&) const;

  void DumpProblem(G4double E, G4double P1, G4double P2, G4double P) const;

  G4FermiPhaseSpaceDecay(const G4FermiPhaseSpaceDecay&);
  const G4FermiPhaseSpaceDecay & operator=(const G4FermiPhaseSpaceDecay &); 
  G4bool operator==(const G4FermiPhaseSpaceDecay&);
  G4bool operator!=(const G4FermiPhaseSpaceDecay&);

  G4Pow* g4pow;
};

inline G4double 
G4FermiPhaseSpaceDecay::PtwoBody(G4double E, G4double P1, G4double P2) const
{
  G4double res = 0.0;
  G4double P = (E+P1+P2)*(E+P1-P2)*(E-P1+P2)*(E-P1-P2)/(4.0*E*E);
  if (P>0.0) { res = std::sqrt(P); }
  else { DumpProblem(E,P1,P2,P); }
  return res;
}

inline std::vector<G4LorentzVector*> * G4FermiPhaseSpaceDecay::
Decay(G4double parent_mass_parameter, const std::vector<G4double>& fragment_masses) const
{
  return KopylovNBodyDecay(parent_mass_parameter,fragment_masses);
}

inline G4double G4FermiPhaseSpaceDecay::BetaKopylov(G4int K) const
{
  G4int N = 3*K - 5;
  G4double xN = G4double(N);
  G4double F;
  //G4double Fmax=std::pow((3.*K-5.)/(3.*K-4.),(3.*K-5.)/2.)*std::sqrt(1-((3.*K-5.)/(3.*K-4.))); 
  // VI variant
  G4double Fmax = std::sqrt(g4pow->powN(xN/(xN + 1),N)/(xN + 1)); 
  G4double chi;
  do {
    chi = G4UniformRand();
    F = std::sqrt(g4pow->powN(chi,N)*(1-chi));      
  } while ( Fmax*G4UniformRand() > F);  
  return chi;
}

inline G4ThreeVector 
G4FermiPhaseSpaceDecay::IsotropicVector(G4double Magnitude) const
  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1.0
{
  G4double CosTheta = 2.0*G4UniformRand() - 1.0;
  G4double SinTheta = std::sqrt((1. - CosTheta)*(1. + CosTheta));
  G4double Phi = CLHEP::twopi*G4UniformRand();
  G4ThreeVector Vector(Magnitude*std::cos(Phi)*SinTheta,
		       Magnitude*std::sin(Phi)*SinTheta,
		       Magnitude*CosTheta);
  return Vector;
}

#endif

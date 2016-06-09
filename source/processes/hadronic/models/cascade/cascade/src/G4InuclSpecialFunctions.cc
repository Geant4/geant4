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
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100914  M. Kelsey -- Migrate to integer A and Z.  Discard pointless
//		verbosity.
// 20120608  M. Kelsey -- Fix variable-name "shadowing" compiler warnings.

#include <cmath>

#include "G4InuclSpecialFunctions.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"


G4double G4InuclSpecialFunctions::getAL(G4int A) {
  return 0.76 + 2.2 / G4cbrt(A);
}

G4double G4InuclSpecialFunctions::csNN(G4double e) {
  G4double snn;

  if (e < 40.0) {
    snn = -1174.8 / (e * e) + 3088.5 / e + 5.3107;
  } else {
    snn = 93074.0 / (e * e) - 11.148 / e + 22.429;
  }

  return snn; 
}

G4double G4InuclSpecialFunctions::csPN(G4double e) {
  G4double spn;

  if (e < 40.0) {
    spn = -5057.4 / (e * e) + 9069.2 / e + 6.9466;
  } else {
    spn = 239380.0 / (e * e) + 1802.0 / e + 27.147;
  }

  return spn; 
}

// calculates the nuclei Fermi energy for 0 - neutron and 1 - proton

G4double G4InuclSpecialFunctions::FermiEnergy(G4int A, G4int Z, G4int ntype) {
  const G4double C = 55.4;
  G4double arg = (ntype==0) ? G4double(A-Z)/A : G4double(Z)/A;

  return C * G4cbrt(arg*arg);	// 2/3 power
}

G4double G4InuclSpecialFunctions::G4cbrt(G4double x) {
  return x==0 ? 0. : (x<0?-1.:1.)*std::exp(std::log(std::fabs(x))/3.);
}

G4double G4InuclSpecialFunctions::inuclRndm() { 
  return G4UniformRand(); 
} 

G4double G4InuclSpecialFunctions::randomGauss(G4double sigma) {
  const G4double eps = 1.0e-6;
  G4double r1 = inuclRndm();
  r1 = r1 > eps ? r1 : eps;
  G4double r2 = inuclRndm();
  r2 = r2 > eps ? r2 : eps;
  r2 = r2 < 1.0 - eps ? r2 : 1.0 - eps; 

  return sigma * std::sin(twopi * r1) * std::sqrt(-2.0 * std::log(r2)); 
} 

G4double G4InuclSpecialFunctions::randomPHI() { 
  return twopi * inuclRndm();
} 

std::pair<G4double, G4double> G4InuclSpecialFunctions::randomCOS_SIN() {
  G4double CT = 1.0 - 2.0 * inuclRndm();

  return std::pair<G4double, G4double>(CT, std::sqrt(1.0 - CT*CT));
}

G4LorentzVector 
G4InuclSpecialFunctions::generateWithFixedTheta(G4double ct, G4double p, 
						G4double mass) {
  G4double phi = randomPHI();
  G4double pt = p * std::sqrt(std::fabs(1.0 - ct * ct));

  static G4ThreeVector pvec;	// Buffers to avoid memory thrashing
  static G4LorentzVector momr;

  pvec.set(pt*std::cos(phi), pt*std::sin(phi), p*ct);
  momr.setVectM(pvec, mass);

  return momr;
}

G4LorentzVector 
G4InuclSpecialFunctions::generateWithRandomAngles(G4double p, G4double mass) {
  std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
  G4double phi = randomPHI();
  G4double pt = p * COS_SIN.second;
  
  static G4ThreeVector pvec;	// Buffers to avoid memory thrashing
  static G4LorentzVector momr;

  pvec.set(pt*std::cos(phi), pt*std::sin(phi), p*COS_SIN.first);
  momr.setVectM(pvec, mass);

  return momr;
}

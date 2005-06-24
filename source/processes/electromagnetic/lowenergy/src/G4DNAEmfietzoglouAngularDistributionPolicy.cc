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
//
// $Id: G4DNAEmfietzoglouAngularDistributionPolicy.cc,v 1.1 2005-06-24 10:07:13 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAEmfietzoglouAngularDistributionPolicy.hh"

#include "Randomize.hh"

G4double                                G4DNAEmfietzoglouAngularDistributionPolicy :: RandomizeCosTheta(G4double k, G4int z) const
{
 //  d sigma_el                sigma_Ruth(K)
 // ------------ (K) ~ -----------------------------  
 //   d Omega           (1 + 2 n(K) - cos(theta))^2
 //
 // We extract cos(theta) distributed as (1 + 2 n(K) - cos(theta))^-2
 //
 // Maximum is for theta=0: 1/(4 n(K)^2) (When n(K) is positive, that is always satisfied within the validity of the process)
 //
 // Phys. Med. Biol. 45 (2000) 3171-3194

 G4double n;
 n=ScreeningFactor(k, z);

 G4double oneOverMax;
 oneOverMax=(4.*n*n);

 G4double cosTheta;
 G4double fCosTheta;

 do
 {
  cosTheta = 2.*G4UniformRand()-1.;
  fCosTheta = (1 + 2.*n - cosTheta);
  fCosTheta = oneOverMax/(fCosTheta*fCosTheta);
 }
 while (fCosTheta < G4UniformRand());

 return cosTheta;
}



G4double                               G4DNAEmfietzoglouAngularDistributionPolicy :: ScreeningFactor(G4double k, G4int z) const
{
 //
 //         alpha_1 + beta_1 ln(K/eV)   constK Z^(2/3)
 // n(T) = -------------------------- -----------------
 //              K/(m_e c^2)            2 + K/(m_e c^2)
 //
 // Where K is the electron non-relativistic kinetic energy
 //
 // n(T) > 0 for T < ~ 400 MeV
 // 
 // Nucl. Instr. Meth. 155 (1978) 145-156

 const G4double alpha_1(1.64);
 const G4double beta_1(-0.0825);
 const G4double constK(1.7E-5);
 
 G4double numerator;
 numerator=(alpha_1+beta_1*std::log(k/eV))*constK*std::pow(static_cast<double>(z), 2./3.);

 k/=electron_mass_c2;

 G4double denominator;
 denominator=k*(2+k);
 
 return numerator/denominator;
}

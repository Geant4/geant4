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
// $Id: G4FinalStateElasticScreenedRutherford.cc,v 1.7 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4FinalStateElasticScreenedRutherford.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateElasticScreenedRutherford::G4FinalStateElasticScreenedRutherford()
{
  lowEnergyLimit = 200 * eV;
  highEnergyLimit = 10 * MeV;

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4FinalStateElasticScreenedRutherford is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateElasticScreenedRutherford::~G4FinalStateElasticScreenedRutherford()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4FinalStateProduct& G4FinalStateElasticScreenedRutherford::GenerateFinalState(const G4Track& track, const G4Step& )
{
  product.Clear();

  G4double k = track.GetDynamicParticle()->GetKineticEnergy();

  const G4int z = 10; 

  G4double cosTheta = RandomizeCosTheta(k, z);
  
  G4double phi = 2. * pi * G4UniformRand();

  G4ThreeVector zVers = track.GetDynamicParticle()->GetMomentumDirection();
  G4ThreeVector xVers = zVers.orthogonal();
  G4ThreeVector yVers = zVers.cross(xVers);

  G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
  G4double yDir = xDir;
  xDir *= std::cos(phi);
  yDir *= std::sin(phi);

  G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));

  product.ModifyPrimaryParticle(zPrimeVers,k);

  return product;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateElasticScreenedRutherford::RandomizeCosTheta(G4double k, G4int z) const
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

 G4double n = ScreeningFactor(k, z);

 G4double oneOverMax = (4. * n*n);

 G4double cosTheta;
 G4double fCosTheta;

 do 
 { 
   cosTheta = 2. * G4UniformRand() - 1.;
   fCosTheta = (1 + 2.*n - cosTheta);
   fCosTheta = oneOverMax / (fCosTheta*fCosTheta);
 }
 while (fCosTheta < G4UniformRand());
 
 return cosTheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateElasticScreenedRutherford::ScreeningFactor(G4double k, G4int z) const
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
  
  const G4double alpha_1 = 1.64;
  const G4double beta_1 = -0.0825;
  const G4double constK = 1.7E-5;
  
  G4double numerator = (alpha_1 + beta_1 * std::log(k/eV)) * constK * std::pow(static_cast<double>(z), 2./3.);
  
  k /= electron_mass_c2;
  
  G4double denominator;
  denominator = k * (2 + k);
  
  G4double result = 0.;
  if (denominator != 0.) 
  {
    result = numerator / denominator;
  }
  else
  {
    G4Exception("G4FinalStateElasticScreenedRutherford::ScreeningFactor - denominator = 0");
  }
  return result;
}

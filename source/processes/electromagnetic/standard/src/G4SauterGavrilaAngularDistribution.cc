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
// $Id: G4SauterGavrilaAngularDistribution.cc 91726 2015-08-03 15:41:36Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4SauterGavrilaAngularDistribution
//
// Author:     Vladimir Ivanchenko using Michel Maire algorithm
//             developed for Geant3
// 
// Creation date: 23 July 2012
//
//
// -------------------------------------------------------------------
//

#include "G4SauterGavrilaAngularDistribution.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

G4SauterGavrilaAngularDistribution::G4SauterGavrilaAngularDistribution()
  : G4VEmAngularDistribution("AngularGenSauterGavrila")
{}    

G4SauterGavrilaAngularDistribution::~G4SauterGavrilaAngularDistribution() 
{}

G4ThreeVector& 
G4SauterGavrilaAngularDistribution::SampleDirection(
       const G4DynamicParticle* dp, G4double, G4int, const G4Material*)
{
  G4double tau = dp->GetKineticEnergy()/electron_mass_c2;
  static const G4double taulimit = 50.0;

  if (tau > taulimit) {
    fLocalDirection = dp->GetMomentumDirection(); 
    // Bugzilla 1120
    // SI on 05/09/2010 as suggested by JG 04/09/10 
  } else {
    // algorithm according Penelope 2008 manual and 
    // F.Sauter Ann. Physik 9, 217(1931); 11, 454(1931). 

    G4double gamma = tau + 1;
    G4double beta  = std::sqrt(tau*(tau + 2))/gamma;
    G4double A     = (1 - beta)/beta;
    G4double Ap2   = A + 2;
    G4double B     = 0.5*beta*gamma*(gamma - 1)*(gamma - 2);
    G4double grej  = 2*(1 + A*B)/A;
    G4double z, g;
    do { 
      G4double q = G4UniformRand();
      z = 2*A*(2*q + Ap2*std::sqrt(q))/(Ap2*Ap2 - 4*q); 
      g = (2 - z)*(1.0/(A + z) + B);

      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(g < G4UniformRand()*grej);
 
    G4double cost = 1 - z;
    G4double sint = std::sqrt(z*(2 - z));
    G4double phi  = CLHEP::twopi*G4UniformRand(); 

    fLocalDirection.set(sint*std::cos(phi), sint*std::sin(phi), cost);
    fLocalDirection.rotateUz(dp->GetMomentumDirection());
  }
  return fLocalDirection;
}

void G4SauterGavrilaAngularDistribution::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Non-polarized photoelectric effect angular generator." << G4endl;
  G4cout << "The Sauter-Gavrila distribution for the K-shell is used."<<G4endl;
  G4cout << "Originally developed by M.Maire for Geant3" 
         << G4endl;
} 

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
// $Id: G4SauterGavrilaAngularDistribution.cc 104021 2017-05-08 07:35:57Z gcosmo $
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
// Modified: 
//   04.05.2017 Marilena Bandieramonte implemented Penelope 2014 algorithm
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

G4ThreeVector&  G4SauterGavrilaAngularDistribution::SampleDirection( 
       const G4DynamicParticle* dp, G4double, G4int, const G4Material*)
{
  static const G4double emin = 1*CLHEP::eV;
  static const G4double emax = 100*CLHEP::MeV;

  G4double energy = std::max(dp->GetKineticEnergy(), emin);
  if (energy > emax) {
    fLocalDirection = dp->GetMomentumDirection();
  } else {
    // Initial algorithm according Penelope 2008 manual and 
    // F.Sauter Ann. Physik 9, 217(1931); 11, 454(1931).
    // Modified according Penelope 2014 manual 
    G4double costheta = 0.0;

    // 1) initialize energy-dependent variables
    // Variable naming according to Eq. (2.24) of Penelope Manual
    // (pag. 44)
    G4double tau = energy/electron_mass_c2;
    G4double gamma = 1.0 + tau;
    G4double beta = std::sqrt(tau*(tau + 2.0))/gamma;
    
    // ac corresponds to "A" of Eq. (2.31)
    //
    G4double ac = (1.0 - beta)/beta;
    G4double a1 = 0.5*beta*gamma*tau*(gamma-2.0);
    G4double a2 = ac + 2.0;
    // gtmax = maximum of the rejection function according to Eq. (2.28), 
    // obtained for tsam=0
    G4double gtmax = 2.0*(a1 + 1.0/ac);
    
    G4double tsam = 0.0;
    G4double gtr  = 0.0;
    
    //2) sampling. Eq. (2.31) of Penelope Manual
    // tsam = 1-std::cos(theta)
    // gtr = rejection function according to Eq. (2.28)
    do{
      G4double rand = G4UniformRand();
      tsam = 2.0*ac * (2.0*rand + a2*std::sqrt(rand)) / (a2*a2 - 4.0*rand);
      gtr = (2.0 - tsam) * (a1 + 1.0/(ac+tsam));
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(G4UniformRand()*gtmax > gtr);

    costheta = 1.0 - tsam;
        
    G4double sint = std::sqrt(tsam*(2.0 - tsam));
    G4double phi  = CLHEP::twopi*G4UniformRand();
        
    fLocalDirection.set(sint*std::cos(phi), sint*std::sin(phi), costheta);
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


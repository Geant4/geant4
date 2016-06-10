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
// $Id: G4BohrFluctuations.cc 91726 2015-08-03 15:41:36Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BohrFluctuations
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 02.04.2003
//
// Modifications:
//
// 23-05-03  Add control on parthalogical cases (V.Ivanchenko)
// 16-10-03 Changed interface to Initialisation (V.Ivanchenko)
//
// Class Description: Sampling of Gaussion fluctuations
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4BohrFluctuations.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4BohrFluctuations::G4BohrFluctuations(const G4String& nam)
 :G4VEmFluctuationModel(nam),
  particle(0),
  minNumberInteractionsBohr(2.0),
  minFraction(0.2),
  xmin(0.2),
  minLoss(0.001*eV)
{
  particleMass   = proton_mass_c2;
  chargeSquare   = 1.0;
  kineticEnergy  = 0.0;
  beta2          = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BohrFluctuations::~G4BohrFluctuations()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BohrFluctuations::InitialiseMe(const G4ParticleDefinition* part)
{
  particle       = part;
  particleMass   = part->GetPDGMass();
  G4double q     = part->GetPDGCharge()/eplus;
  chargeSquare   = q*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4BohrFluctuations::SampleFluctuations(const G4MaterialCutsCouple* couple,
                                       const G4DynamicParticle* dp,
                                       G4double tmax,
                                       G4double length,
                                       G4double meanLoss)
{
  if(meanLoss <= minLoss) { return meanLoss; }
  const G4Material* material = couple->GetMaterial();
  G4double siga = Dispersion(material,dp,tmax,length);
  G4double loss = meanLoss;

  G4double navr = meanLoss*meanLoss/siga;
  //G4cout << "### meanLoss= " << meanLoss << "  navr= " << navr << G4endl;
  if (navr >= minNumberInteractionsBohr) {
 
    // Increase fluctuations for big fractional energy loss
    if ( meanLoss > minFraction*kineticEnergy ) {
      G4double gam = (kineticEnergy - meanLoss)/particleMass + 1.0;
      G4double b2  = 1.0 - 1.0/(gam*gam);
      if(b2 < xmin*beta2) b2 = xmin*beta2;
      G4double x   = b2/beta2;
      G4double x3  = 1.0/(x*x*x);
      siga *= 0.25*(1.0 + x)*(x3 + (1.0/b2 - 0.5)/(1.0/beta2 - 0.5) );
    }
    siga = sqrt(siga);
    G4double twomeanLoss = meanLoss + meanLoss;
    //G4cout << "siga= " << siga << " 2edp= " << twomeanLoss <<G4endl;

    if(twomeanLoss < siga) {
      G4double x;
      do {
        loss = twomeanLoss*G4UniformRand();
        x = (loss - meanLoss)/siga;
        // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
      } while (1.0 - 0.5*x*x < G4UniformRand());
    } else {
      do {
        loss = G4RandGauss::shoot(meanLoss,siga);
        // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
      } while (0.0 > loss || loss > twomeanLoss);
    }

  // Poisson fluctuations
  } else {
    G4double n    = (G4double)(G4Poisson(navr));
    loss = meanLoss*n/navr;
  }
  //G4cout << "loss= " << loss << G4endl;

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BohrFluctuations::Dispersion(const G4Material* material,
                                        const G4DynamicParticle* dp,
                                        G4double tmax,
                                        G4double length)
{
  if(!particle) { InitialiseMe(dp->GetDefinition()); }

  G4double electronDensity = material->GetElectronDensity();
  kineticEnergy = dp->GetKineticEnergy();
  G4double etot = kineticEnergy + particleMass;
  beta2 = kineticEnergy*(kineticEnergy + 2.0*particleMass)/(etot*etot);
  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                 * electronDensity * chargeSquare;

  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



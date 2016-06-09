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
// $Id: G4BohrFluctuations.cc,v 1.2 2004/12/01 19:37:14 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
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
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4BohrFluctuations::G4BohrFluctuations(const G4String& nam)
 :G4VEmFluctuationModel(nam),
  particle(0),
  minNumberInteractionsBohr(10.0),
  minFraction(0.2),
  xmin(0.2),
  minLoss(0.001*eV)
{}

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

G4double G4BohrFluctuations::SampleFluctuations(const G4Material* material,
                                                const G4DynamicParticle* dp,
                                                      G4double& tmax,
                                                      G4double& length,
                                                      G4double& meanLoss)
{
  if(meanLoss <= minLoss) return meanLoss;
  G4double siga = Dispersion(material,dp,tmax,length);
  G4double loss = meanLoss;

  G4double navr = minNumberInteractionsBohr;

  // Gaussian fluctuation
  G4bool gauss = true;
  if (meanLoss < minNumberInteractionsBohr*tmax) {
    navr = meanLoss*meanLoss/siga;
    if (navr < minNumberInteractionsBohr) gauss = false;
  }
//  G4cout << "### meanLoss= " << meanLoss << "  navr= " << navr << " sig= " << sqrt(siga) << G4endl;

  if(gauss) {
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

    if(twomeanLoss < siga) {
      G4double x;
      do {
        loss = twomeanLoss*G4UniformRand();
	x = (loss - meanLoss)/siga;
      } while (1.0 - 0.5*x*x < G4UniformRand());
    } else {
      do {
        loss = G4RandGauss::shoot(meanLoss,siga);
      } while (0.0 > loss || loss > twomeanLoss);
    }

  // Poisson fluctuations
  } else {
    G4double n    = (G4double)(G4Poisson(navr));
    loss = meanLoss*n/navr;
  }

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4DynamicParticleFluctuation
//
// Author:        V. Ivanchenko 
// 
// Creation date: 23.08.2024
//
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4DynamicParticleFluctuation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4Log.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DynamicParticleFluctuation::G4DynamicParticleFluctuation(const G4String& nam)
 : G4UniversalFluctuation(nam)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DynamicParticleFluctuation::InitialiseLocal(const G4DynamicParticle* part)
{
  particleMass = part->GetMass();
  const G4double q = part->GetCharge()/CLHEP::eplus;

  // Derived quantities
  m_Inv_particleMass = 1.0 / particleMass;
  m_massrate = CLHEP::electron_mass_c2 * m_Inv_particleMass;
  chargeSquare = q*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DynamicParticleFluctuation::SampleFluctuations(
					 const G4MaterialCutsCouple* couple,
                                         const G4DynamicParticle* dp,
                                         const G4double tcut,
                                         const G4double tmax,
                                         const G4double length,
                                         const G4double averageLoss)
{
  // Calculate actual loss from the mean loss.
  // The model used to get the fluctuations is essentially the same
  // as in Glandz in Geant3 (Cern program library W5013, phys332).
  // L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

  // shortcut for very small loss or from a step nearly equal to the range
  // (out of validity of the model)
  //
  if (averageLoss < minLoss) { return averageLoss; }
  meanLoss = averageLoss;
  const G4double tkin = dp->GetKineticEnergy();
  //G4cout<< "Emean= "<< meanLoss<< " tmax= "<< tmax<< " L= "<<length<<G4endl;

  CLHEP::HepRandomEngine* rndmEngineF = G4Random::getTheEngine();

  InitialiseLocal(dp);             
  const G4double gam   = tkin * m_Inv_particleMass + 1.0;
  const G4double gam2  = gam*gam;
  const G4double beta  = dp->GetBeta(); 
  const G4double beta2 = beta*beta;

  G4double loss(0.), siga(0.);

  const G4Material* material = couple->GetMaterial();
  
  // Gaussian regime
  // for heavy particles only and conditions
  // for Gauusian fluct. has been changed 
  //
  if (particleMass > CLHEP::electron_mass_c2 &&
      meanLoss >= minNumberInteractionsBohr*tcut && tmax <= 2.*tcut) {

    siga = std::sqrt((tmax/beta2 - 0.5*tcut)*CLHEP::twopi_mc2_rcl2* 
                      length*chargeSquare*material->GetElectronDensity());
    const G4double sn = meanLoss/siga;
  
    // thick target case 
    if (sn >= 2.0) {

      const G4double twomeanLoss = meanLoss + meanLoss;
      do {
	loss = G4RandGauss::shoot(rndmEngineF, meanLoss, siga);
	// Loop checking, 03-Aug-2015, Vladimir Ivanchenko
      } while  (0.0 > loss || twomeanLoss < loss);

      // Gamma distribution
    } else {

      const G4double neff = sn*sn;
      loss = meanLoss*G4RandGamma::shoot(rndmEngineF, neff, 1.0)/neff;
    }
    //G4cout << "Gauss: " << loss << G4endl;
    return loss;
  }

  auto ioni = material->GetIonisation();
  e0 = ioni->GetEnergy0fluct();

  // very small step or low-density material
  if(tcut <= e0) { return meanLoss; }

  ipotFluct = ioni->GetMeanExcitationEnergy();
  ipotLogFluct = ioni->GetLogMeanExcEnergy();

  // width correction for small cuts
  const G4double scaling = std::min(1.+0.5*CLHEP::keV/tcut, 1.50);
  meanLoss /= scaling;

  w2 = (tcut > ipotFluct) ? 
    G4Log(2.*CLHEP::electron_mass_c2*beta2*gam2) - beta2 : 0.0;
  return SampleGlandz(rndmEngineF, material, tcut)*scaling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4DynamicParticleFluctuation::Dispersion(
                          const G4Material* material,
                          const G4DynamicParticle* dp,
                          const G4double tcut,
                          const G4double tmax,
                          const G4double length)
{
  InitialiseLocal(dp);
  const G4double beta = dp->GetBeta();
  return (tmax/(beta*beta) - 0.5*tcut) * CLHEP::twopi_mc2_rcl2 * length
    * material->GetElectronDensity() * chargeSquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

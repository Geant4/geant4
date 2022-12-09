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
// File name:     G4UniversalFluctuation
//
// Author:        V. Ivanchenko for Laszlo Urban
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UniversalFluctuation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Log.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UniversalFluctuation::G4UniversalFluctuation(const G4String& nam)
 :G4VEmFluctuationModel(nam),
  minLoss(10.*CLHEP::eV)
{
  rndmarray = new G4double[sizearray];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UniversalFluctuation::~G4UniversalFluctuation()
{
  delete [] rndmarray;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UniversalFluctuation::InitialiseMe(const G4ParticleDefinition* part)
{
  particle = part;
  particleMass = part->GetPDGMass();
  const G4double q = part->GetPDGCharge()/CLHEP::eplus;

  // Derived quantities
  m_Inv_particleMass = 1.0 / particleMass;
  m_massrate = CLHEP::electron_mass_c2 * m_Inv_particleMass;
  chargeSquare = q*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4UniversalFluctuation::SampleFluctuations(const G4MaterialCutsCouple* couple,
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
  const G4double tkin  = dp->GetKineticEnergy();
  //G4cout<< "Emean= "<< meanLoss<< " tmax= "<< tmax<< " L= "<<length<<G4endl;

  if(dp->GetDefinition() != particle) { InitialiseMe(dp->GetDefinition()); }

  CLHEP::HepRandomEngine* rndmEngineF = G4Random::getTheEngine();
             
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
    G4Log(2.*CLHEP::electron_mass_c2*beta2*gam2)-beta2 : 0.0;
  return SampleGlandz(rndmEngineF, material, tcut)*scaling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4UniversalFluctuation::SampleGlandz(CLHEP::HepRandomEngine* rndmEngineF,
                                     const G4Material*,
                                     const G4double tcut)
{
  G4double a1(0.0), a3(0.0);
  G4double loss = 0.0;
  G4double e1 = ipotFluct;

  if(tcut > e1) {
    a1 = meanLoss*(1.-rate)/e1;
    if(a1 < a0) {
      const G4double fwnow = 0.1+(fw-0.1)*std::sqrt(a1/a0);
      a1 /= fwnow;
      e1 *= fwnow;
    } else {
      a1 /= fw;
      e1 *= fw;
    }   
  }

  const G4double w1 = tcut/e0;
  a3 = rate*meanLoss*(tcut - e0)/(e0*tcut*G4Log(w1));
  if(a1 <= 0.) { a3 /= rate; }
  
  //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont  
  G4double emean = 0.;
  G4double sig2e = 0.;

  // excitation of type 1
  if(a1 > 0.0) { AddExcitation(rndmEngineF, a1, e1, emean, loss, sig2e); }

  if(sig2e > 0.0) { SampleGauss(rndmEngineF, emean, sig2e, loss); }

  // ionisation 
  if(a3 > 0.) {
    emean = 0.;
    sig2e = 0.;
    G4double p3 = a3;
    G4double alfa = 1.;
    if(a3 > nmaxCont) {
      alfa = w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
      const G4double alfa1  = alfa*G4Log(alfa)/(alfa-1.);
      const G4double namean = a3*w1*(alfa-1.)/((w1-1.)*alfa);
      emean += namean*e0*alfa1;
      sig2e += e0*e0*namean*(alfa-alfa1*alfa1);
      p3 = a3 - namean;
    }

    const G4double w3 = alfa*e0;
    if(tcut > w3) {
      const G4double w = (tcut-w3)/tcut;
      const G4int nnb = (G4int)G4Poisson(p3);
      if(nnb > 0) {
        if(nnb > sizearray) {
          sizearray = nnb;
          delete [] rndmarray;
          rndmarray = new G4double[nnb];
        }
        rndmEngineF->flatArray(nnb, rndmarray);
        for (G4int k=0; k<nnb; ++k) { loss += w3/(1.-w*rndmarray[k]); }
      }
    }
    if(sig2e > 0.0) { SampleGauss(rndmEngineF, emean, sig2e, loss); }
  }
  //G4cout << "### loss=" << loss << G4endl;
  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4UniversalFluctuation::Dispersion(
                          const G4Material* material,
                          const G4DynamicParticle* dp,
                          const G4double tcut,
                          const G4double tmax,
                          const G4double length)
{
  if(dp->GetDefinition() != particle) { InitialiseMe(dp->GetDefinition()); }
  const G4double beta = dp->GetBeta();
  return (tmax/(beta*beta) - 0.5*tcut) * CLHEP::twopi_mc2_rcl2 * length
    * material->GetElectronDensity() * chargeSquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4UniversalFluctuation::SetParticleAndCharge(const G4ParticleDefinition* part,
                                             G4double q2)
{
  if(part != particle) {
    particle = part;
    particleMass = part->GetPDGMass();

    // Derived quantities
    m_Inv_particleMass = 1.0 / particleMass;
    m_massrate = CLHEP::electron_mass_c2 * m_Inv_particleMass;
  }
  chargeSquare = q2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

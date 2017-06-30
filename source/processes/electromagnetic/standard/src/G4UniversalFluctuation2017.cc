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
// $Id: G4UniversalFluctuation2017.cc 104232 2017-05-19 14:56:13Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4UniversalFluctuation2017
//
// Author:        V. Ivanchenko for Laszlo Urban
// 
// Creation date: 24.05.2017
//
// Modifications: 
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UniversalFluctuation2017.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4UniversalFluctuation2017::G4UniversalFluctuation2017(const G4String& nam)
 :G4VEmFluctuationModel(nam),
  particle(nullptr),
  minNumberInteractionsBohr(10.0),
  minLoss(10.*eV),
  nmaxCont(16.),
  rate(0.511),
  fw(3.)
{
  lastMaterial = nullptr;
  fw2 = fw*fw;
  fw3 = (fw-1.)/(fw2*fw2);
  particleMass = chargeSquare = ipotFluct = electronDensity = f1Fluct = f2Fluct 
    = e1Fluct = e2Fluct = e1LogFluct = e2LogFluct = ipotLogFluct = e0 = esmall 
    = e1 = e2 = 0.0;
  m_Inv_particleMass = m_massrate = DBL_MAX;
  sizearray = 30;
  rndmarray = new G4double[30];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UniversalFluctuation2017::~G4UniversalFluctuation2017()
{
  delete [] rndmarray;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UniversalFluctuation2017::InitialiseMe(const G4ParticleDefinition* part)
{
  particle       = part;
  particleMass   = part->GetPDGMass();
  G4double q     = part->GetPDGCharge()/eplus;

  // Derived quantities
  m_Inv_particleMass = 1.0 / particleMass;
  m_massrate = electron_mass_c2 * m_Inv_particleMass ;
  chargeSquare   = q*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4UniversalFluctuation2017::SampleFluctuations(const G4MaterialCutsCouple* couple,
                                           const G4DynamicParticle* dp,
                                           G4double tmax,
                                           G4double length,
                                           G4double averageLoss)
{
  // Calculate actual loss from the mean loss.
  // The model used to get the fluctuations is essentially the same
  // as in Glandz in Geant3 (Cern program library W5013, phys332).
  // L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

  // shortcut for very small loss or from a step nearly equal to the range
  // (out of validity of the model)
  //
  G4double meanLoss = averageLoss;
  G4double tkin  = dp->GetKineticEnergy();
  //G4cout<< "Emean= "<< meanLoss<< " tmax= "<< tmax<< " L= "<<length<<G4endl;
  if (meanLoss < minLoss) { return meanLoss; }

  if(dp->GetDefinition() != particle) { InitialiseMe(dp->GetDefinition()); }

  CLHEP::HepRandomEngine* rndmEngineF = G4Random::getTheEngine();
  
  G4double tau   = tkin * m_Inv_particleMass;            
  G4double gam   = tau + 1.0;
  G4double gam2  = gam*gam;
  G4double beta2 = tau*(tau + 2.0)/gam2;

  G4double loss(0.), siga(0.);

  const G4Material* material = couple->GetMaterial();
  
  // Gaussian regime
  // for heavy particles only and conditions
  // for Gauusian fluct. has been changed 
  //
  if ((particleMass > electron_mass_c2) &&
      (meanLoss >= minNumberInteractionsBohr*tmax))
  {
    G4double tmaxkine = 2.*electron_mass_c2*beta2*gam2/
                        (1.+m_massrate*(2.*gam+m_massrate)) ;
    if (tmaxkine <= 2.*tmax)   
    {
      electronDensity = material->GetElectronDensity();
      siga = sqrt((1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                  * electronDensity * chargeSquare);

      G4double sn = meanLoss/siga;
  
      // thick target case 
      if (sn >= 2.0) {

        G4double twomeanLoss = meanLoss + meanLoss;
        do {
          loss = G4RandGauss::shoot(rndmEngineF,meanLoss,siga);
          // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
        } while  (0.0 > loss || twomeanLoss < loss);

        // Gamma distribution
      } else {

        G4double neff = sn*sn;
        loss = meanLoss*G4RandGamma::shoot(rndmEngineF,neff,1.0)/neff;
      }
      //G4cout << "Gauss: " << loss << G4endl;
      return loss;
    }
  }

  // Glandz regime : initialisation
  //
  if (material != lastMaterial) {
    f1Fluct      = material->GetIonisation()->GetF1fluct();
    f2Fluct      = material->GetIonisation()->GetF2fluct();
    e1Fluct      = material->GetIonisation()->GetEnergy1fluct();
    e2Fluct      = material->GetIonisation()->GetEnergy2fluct();
    e1LogFluct   = material->GetIonisation()->GetLogEnergy1fluct();
    e2LogFluct   = material->GetIonisation()->GetLogEnergy2fluct();
    ipotFluct    = material->GetIonisation()->GetMeanExcitationEnergy();
    ipotLogFluct = material->GetIonisation()->GetLogMeanExcEnergy();
    e0 = material->GetIonisation()->GetEnergy0fluct();
    esmall = 0.5*sqrt(e0*ipotFluct);  
    lastMaterial = material;   
  }

  // very small step or low-density material
  if(tmax <= e0) { return meanLoss; }

  G4double a1(0.0), a2(0.0), a3(0.0);
    
  loss = 0.0;

  G4double rnow = rate + 0.489/(1.+0.050*meanLoss/ipotFluct);
  e1 = e1Fluct;
  e2 = e2Fluct;

  if(tmax > ipotFluct) {
    G4double w2 = G4Log(2.*electron_mass_c2*beta2*gam2)-beta2;

    if(w2 > ipotLogFluct)  {
      if(w2 > e2LogFluct) {
	G4double C = meanLoss*(1.-rnow)/(w2-ipotLogFluct);
	a1 = C*f1Fluct*(w2-e1LogFluct)/e1Fluct;
	a2 = C*f2Fluct*(w2-e2LogFluct)/e2Fluct;
      } else {
	a1 = meanLoss*(1.-rnow)/e1;
      }
      if(a1 < fw2) { 
        G4double fwnow = fw-fw3*(fw2-a1)*(fw2-a1);
        a1 /= fwnow;
        e1 *= fwnow;
      } else {
        a1 /= fw;
        e1 = fw*e1Fluct;
      }
      if(a2 < fw2) { 
        G4double fwnow = fw-fw3*(fw2-a2)*(fw2-a2);
        a2 /= fwnow;
        e2 *= fwnow;
      } else {
        a2 /= fw;
        e2 = fw*e2Fluct;
      }
    }   
  }

  G4double w1 = tmax/e0;
  if(tmax > e0) {
    a3 = rnow*meanLoss*(tmax-e0)/(e0*tmax*G4Log(w1));
    if(a1+a2 <= 0.) { 
      a3 /= rnow; 
    }
  }
  //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont  
  G4double emean = 0.;
  G4double sig2e = 0.;

  // excitation of type 1
  if(a1 > 0.0) { AddExcitation(rndmEngineF, a1, e1, emean, loss, sig2e); }

  // excitation of type 2
  if(a2 > 0.0) { AddExcitation(rndmEngineF, a2, e2, emean, loss, sig2e); }

  if(emean > 0.0) { SampleGauss(rndmEngineF, emean, sig2e, loss); }

  // ionisation 
  if(a3 > 0.) {
    emean = 0.;
    sig2e = 0.;
    G4double p3 = a3;
    G4double alfa = 1.;
    if(a3 > nmaxCont)
      {
        alfa            = w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
        G4double alfa1  = alfa*G4Log(alfa)/(alfa-1.);
        G4double namean = a3*w1*(alfa-1.)/((w1-1.)*alfa);
        emean          += namean*e0*alfa1;
        sig2e          += e0*e0*namean*(alfa-alfa1*alfa1);
        p3              = a3-namean;
      }

    G4double w2 = alfa*e0;
    if(tmax > w2) {
      G4double w  = (tmax-w2)/tmax;
      G4int nnb = G4Poisson(p3);
      if(nnb > 0) {
        if(nnb > sizearray) {
          sizearray = nnb;
          delete [] rndmarray;
          rndmarray = new G4double[nnb];
        }
        rndmEngineF->flatArray(nnb, rndmarray);
        for (G4int k=0; k<nnb; ++k) { loss += w2/(1.-w*rndmarray[k]); }
      }
    }
    if(emean > 0.0) { SampleGauss(rndmEngineF, emean, sig2e, loss); }
  }

  return loss;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4UniversalFluctuation2017::Dispersion(
                          const G4Material* material,
                          const G4DynamicParticle* dp,
                                G4double tmax,
                                G4double length)
{
  if(dp->GetDefinition() != particle) { InitialiseMe(dp->GetDefinition()); }

  electronDensity = material->GetElectronDensity();

  G4double gam   = (dp->GetKineticEnergy())*m_Inv_particleMass + 1.0;
  G4double beta2 = 1.0 - 1.0/(gam*gam);

  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                 * electronDensity * chargeSquare;

  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4UniversalFluctuation2017::SetParticleAndCharge(const G4ParticleDefinition* part,
                                             G4double q2)
{
  if(part != particle) {
    particle       = part;
    particleMass   = part->GetPDGMass();

    // Derived quantities
    if( particleMass != 0.0 ){
      m_Inv_particleMass = 1.0 / particleMass;
      m_massrate = electron_mass_c2 * m_Inv_particleMass ;
    }else{
      m_Inv_particleMass = DBL_MAX;
      m_massrate = DBL_MAX;
    }
  }
  chargeSquare = q2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

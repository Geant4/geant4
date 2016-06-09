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
// $Id: G4UniversalFluctuation.cc,v 1.5 2005/06/28 09:13:23 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01-patch-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4UniversalFluctuation
//
// Author:        Vladimir Ivanchenko 
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
// 28-12-02 add method Dispersion (V.Ivanchenko)
// 07-02-03 change signature (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 16-10-03 Changed interface to Initialisation (V.Ivanchenko)
// 07-11-03 Fix problem of rounding of double in G4UniversalFluctuations
// 06-02-04 Add control on big sigma > 2*meanLoss (V.Ivanchenko)
// 26-04-04 Comment out the case of very small step (V.Ivanchenko)
// 07-02-05 define problim = 5.e-3 (mma)
// 03-05-05 conditions of Gaussian fluctuation changed (bugfix)
//          + smearing for very small loss (L.Urban)
//          

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UniversalFluctuation.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4UniversalFluctuation::G4UniversalFluctuation(const G4String& nam)
 :G4VEmFluctuationModel(nam),
  particle(0),
  minNumberInteractionsBohr(10.0),
  theBohrBeta2(50.0*keV/proton_mass_c2),
  minLoss(10.*eV),
  problim(5.e-3),  
  alim(10.),
  nmaxCont1(4.),
  nmaxCont2(16.)
{
  sumalim = -log(problim);
  lastMaterial = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UniversalFluctuation::~G4UniversalFluctuation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UniversalFluctuation::InitialiseMe(const G4ParticleDefinition* part)
{
  particle       = part;
  particleMass   = part->GetPDGMass();
  G4double q     = part->GetPDGCharge()/eplus;
  chargeSquare   = q*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UniversalFluctuation::SampleFluctuations(const G4Material* material,
                                                const G4DynamicParticle* dp,
				                      G4double& tmax,
					              G4double& length,
                                                      G4double& meanLoss)
{
// Calculate actual loss from the mean loss.
// The model used to get the fluctuations is essentially the same
// as in Glandz in Geant3 (Cern program library W5013, phys332).
// L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

  // shortcut for very very small loss (out of validity of the model)
  //
  if (meanLoss < minLoss) return meanLoss;

  if(!particle) InitialiseMe(dp->GetDefinition());

  G4double tau   = dp->GetKineticEnergy()/particleMass;
  G4double gam   = tau + 1.0;
  G4double gam2  = gam*gam;
  G4double beta2 = tau*(tau + 2.0)/gam2;

  G4double loss(0.), siga(0.);
  
  // Gaussian regime
  // for heavy particles only and conditions
  // for Gauusian fluct. has been changed 
  //
  if ((particleMass > electron_mass_c2) &&
      (meanLoss >= minNumberInteractionsBohr*tmax))
  {
    G4double massrate = electron_mass_c2/particleMass ;
    G4double tmaxkine = 2.*electron_mass_c2*beta2*gam2/
                        (1.+massrate*(2.*gam+massrate)) ;
    if (tmaxkine <= 2.*tmax)   
    {
      electronDensity = material->GetElectronDensity();
      siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                                * electronDensity * chargeSquare;
      siga = sqrt(siga);
      G4double twomeanLoss = meanLoss + meanLoss;
      if (twomeanLoss < siga) {
        G4double x;
        do {
          loss = twomeanLoss*G4UniformRand();
       	  x = (loss - meanLoss)/siga;
        } while (1.0 - 0.5*x*x < G4UniformRand());
      } else {
        do {
          loss = G4RandGauss::shoot(meanLoss,siga);
        } while (loss < 0. || loss > twomeanLoss);
      }
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
    rateFluct    = material->GetIonisation()->GetRateionexcfluct();
    ipotFluct    = material->GetIonisation()->GetMeanExcitationEnergy();
    ipotLogFluct = material->GetIonisation()->GetLogMeanExcEnergy();
    lastMaterial = material;
  }

  G4double a1 = 0. , a2 = 0., a3 = 0. ;
  G4double p1,p2,p3;
  G4double rate = rateFluct ;

  G4double w1 = tmax/ipotFluct;
  G4double w2 = log(2.*electron_mass_c2*beta2*gam2)-beta2;

  if(w2 > ipotLogFluct)
  {
    G4double C = meanLoss*(1.-rateFluct)/(w2-ipotLogFluct);
    a1 = C*f1Fluct*(w2-e1LogFluct)/e1Fluct;
    a2 = C*f2Fluct*(w2-e2LogFluct)/e2Fluct;
    if(a2 < 0.)
    {
      a1 = 0. ;
      a2 = 0. ;
      rate = 1. ;  
    }
  }
  else
  {
    rate = 1. ;
  }

  if(tmax > ipotFluct) 
    a3 = rate*meanLoss*(tmax-ipotFluct)/(ipotFluct*tmax*log(w1));

  G4double suma = a1+a2+a3;
  
  // Glandz regime
  //
  if (suma > sumalim)
  {
    p1 = 0., p2 = 0 ;
    if((a1+a2) > 0.)
    {
      // excitation type 1
      if (a1>alim) {
        siga=sqrt(a1) ;
        p1 = max(0.,G4RandGauss::shoot(a1,siga)+0.5);
      } else {
        p1 = G4double(G4Poisson(a1));
      }
    
      // excitation type 2
      if (a2>alim) {
        siga=sqrt(a2) ;
        p2 = max(0.,G4RandGauss::shoot(a2,siga)+0.5);
      } else {
        p2 = G4double(G4Poisson(a2));
      }
    
      loss = p1*e1Fluct+p2*e2Fluct;
 
      // smearing to avoid unphysical peaks
      if (p2 > 0.)
        loss += (1.-2.*G4UniformRand())*e2Fluct;
      else if (loss>0.)
        loss += (1.-2.*G4UniformRand())*e1Fluct;   
      if (loss < 0.) loss = 0.0;
    }

    // ionisation
    if (a3 > 0.) {
      if (a3>alim) {
        siga=sqrt(a3) ;
        p3 = max(0.,G4RandGauss::shoot(a3,siga)+0.5);
      } else {
        p3 = G4double(G4Poisson(a3));
      }
      G4double lossc = 0.;
      if (p3 > 0) {
        G4double na = 0.; 
        G4double alfa = 1.;
        if (p3 > nmaxCont2) {
          G4double rfac   = p3/(nmaxCont2+p3);
          G4double namean = p3*rfac;
          G4double sa     = nmaxCont1*rfac;
          na              = G4RandGauss::shoot(namean,sa);
          if (na > 0.) {
            alfa   = w1*(nmaxCont2+p3)/(w1*nmaxCont2+p3);
            G4double alfa1  = alfa*log(alfa)/(alfa-1.);
            G4double ea     = na*ipotFluct*alfa1;
            G4double sea    = ipotFluct*sqrt(na*(alfa-alfa1*alfa1));
            lossc += G4RandGauss::shoot(ea,sea);
          }
        }

        if (p3 > na) {
          w2 = alfa*ipotFluct;
          G4double w  = (tmax-w2)/tmax;
          G4int    nb = G4int(p3-na);
          for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
        }
      }        
      loss += lossc;  
    }
    return loss;
  }
  
  // suma < sumalim;  very small energy loss;  
  //
  G4double e0 = material->GetIonisation()->GetEnergy0fluct();

  if(tmax <= e0) a3 = 0.0;
  else a3 = meanLoss*(tmax-e0)/(tmax*e0*log(tmax/e0));

  if (a3 > alim)
  {
    siga=sqrt(a3);
    p3 = max(0.,G4RandGauss::shoot(a3,siga)+0.5);
  } else {
    p3 = G4double(G4Poisson(a3));
  }
  if (p3 > 0.) {
    G4double w = (tmax-e0)/tmax;
    G4double corrfac = 1.;
    if (p3 > nmaxCont2) {
      corrfac = p3/nmaxCont2;
      p3 = nmaxCont2;
    } 
    G4int ip3 = (G4int)p3;
    for (G4int i=0; i<ip3; i++) loss += 1./(1.-w*G4UniformRand());
    loss *= e0*corrfac;
    // smearing for losses near to e0
    if(p3 <= 2.)
    loss += e0*(1.-2.*G4UniformRand()) ;
   }
    
   return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4UniversalFluctuation::Dispersion(
                          const G4Material* material,
                          const G4DynamicParticle* dp,
 				G4double& tmax,
			        G4double& length)
{
  if(!particle) InitialiseMe(dp->GetDefinition());

  electronDensity = material->GetElectronDensity();

  G4double gam   = (dp->GetKineticEnergy())/particleMass + 1.0;
  G4double beta2 = 1.0 - 1.0/(gam*gam);

  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                 * electronDensity * chargeSquare;

  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

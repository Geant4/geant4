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
// Class Description: 
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
#include "G4UniversalFluctuation.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4UniversalFluctuation::G4UniversalFluctuation()
 :G4VEmFluctuationModel(),
  minNumberInteractionsBohr(10.0),
  theBohrBeta2(50.0*keV/proton_mass_c2),
  minLoss(0.000001*eV),
  problim(0.01),
  alim(10.),
  nmaxCont1(4),
  nmaxCont2(16)
{
  lastMaterial = 0;
  sumalim = -log(problim);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4UniversalFluctuation::~G4UniversalFluctuation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UniversalFluctuation::Initialise(const G4ParticleDefinition* part)
{
  particle = part; 
  particleMass = part->GetPDGMass();
  G4double q = part->GetPDGCharge()/eplus;
  chargeSquare = q*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UniversalFluctuation::SampleFluctuations(const G4Material* material, 
                                                    const G4DynamicParticle* dp,
				                          G4double& tmax,
						          G4double& length,
                                                          G4double meanLoss)
{
//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is essentially the same 
// as in Glandz in Geant3.

  // shortcut for very very small loss 
  if(meanLoss < minLoss) return;

  G4double preStepKinEnergy = dp->GetKineticEnergy(); 

  if(dp->GetDefinition() != particle) {
    particleMass = dp->GetMass();
    G4double q = dp->GetCharge();
    chargeSquare = q*q;
  }

  // data members for a given material

  if(material != lastMaterial) {
    ipotFluct = material->GetIonisation()->GetMeanExcitationEnergy();
    electronDensity = material->GetElectronDensity();
    //    zeff = electronDensity/(material->GetTotNbOfAtomsPerVolume());
  }

  // get particle data

  G4double tau   = preStepKinEnergy/particleMass; 
  G4double tau1  = tau + 1.0; 
  G4double tau2  = tau * (tau+2.0);
  G4double beta2 = tau2/(tau1*tau1);


  // Validity range for delta electron cross section
  G4double loss, siga;

  // Gaussian fluctuation 
  if(meanLoss >= minNumberInteractionsBohr*tmax || tmax <= ipotFluct*minNumberInteractionsBohr)
  {
    siga = tmax * (1.0-0.5*beta2) * length * twopi_mc2_rcl2 
         * electronDensity / beta2;
    siga = sqrt(siga * chargeSquare);
/*
    // High velocity or negatively charged particle
    if( beta2 > 3.0*theBohrBeta2*zeff || charge < 0.0) {
      siga = sqrt( siga * chargeSquare ) ;

    // Low velocity - additional ion charge fluctuations according to
    // Q.Yang et al., NIM B61(1991)149-155.
    } else {
      G4double chu = theIonChuFluctuationModel->TheValue(particle, material);
      G4double yang = theIonYangFluctuationModel->TheValue(particle, material);
      siga = sqrt( siga * (chargeSquare * chu + yang)) ;
    }
*/
    do {
     loss = G4RandGauss::shoot(meanLoss,siga);
    } while (loss < 0.);

    meanLoss = loss;
    if(lastMaterial != material) lastMaterial = material;
    return;
  }

  // Non Gaussian fluctuation 

  if(material != lastMaterial) {
    zeff = electronDensity/(material->GetTotNbOfAtomsPerVolume());
    f1Fluct     = material->GetIonisation()->GetF1fluct();
    f2Fluct     = material->GetIonisation()->GetF2fluct();
    e1Fluct     = material->GetIonisation()->GetEnergy1fluct();
    e2Fluct     = material->GetIonisation()->GetEnergy2fluct();
    e1LogFluct  = material->GetIonisation()->GetLogEnergy1fluct();
    e2LogFluct  = material->GetIonisation()->GetLogEnergy2fluct();
    rateFluct   = material->GetIonisation()->GetRateionexcfluct();
    ipotLogFluct= material->GetIonisation()->GetLogMeanExcEnergy();
    lastMaterial = material;
  }

  G4double suma,w1,w2,C,e0,lossc,w;
  G4double a1,a2,a3;
  G4int p1,p2,p3;
  G4int nb;
  G4double corrfac, na,alfa,rfac,namean,sa,alfa1,ea,sea;
  G4double dp3;

  w1 = tmax/ipotFluct;
  w2 = log(2.*electron_mass_c2*tau2);

  C = meanLoss*(1.-rateFluct)/(w2-ipotLogFluct-beta2);

  a1 = C*f1Fluct*(w2-e1LogFluct-beta2)/e1Fluct;
  a2 = C*f2Fluct*(w2-e2LogFluct-beta2)/e2Fluct;
  a3 = rateFluct*meanLoss*(tmax-ipotFluct)/(ipotFluct*tmax*log(w1));
  if(a1 < 0.) a1 = 0.;
  if(a2 < 0.) a2 = 0.;
  if(a3 < 0.) a3 = 0.;

  suma = a1+a2+a3;

  loss = 0. ;

  if(suma < sumalim)             // very small Step
    {
      e0 = material->GetIonisation()->GetEnergy0fluct();

      if(tmax == ipotFluct)
      {
        a3 = meanLoss/e0;

        if(a3>alim)
        {
          siga=sqrt(a3) ;
          p3 = G4std::max(0,int(G4RandGauss::shoot(a3,siga)+0.5));
        }
        else
          p3 = G4Poisson(a3);

        loss = p3*e0 ;

        if(p3 > 0)
          loss += (1.-2.*G4UniformRand())*e0 ;

      }
      else
      {
        tmax = tmax-ipotFluct+e0 ;
        a3 = meanLoss*(tmax-e0)/(tmax*e0*log(tmax/e0));

        if(a3>alim)
        {
          siga=sqrt(a3) ;
          p3 = G4std::max(0,int(G4RandGauss::shoot(a3,siga)+0.5));
        }
        else
          p3 = G4Poisson(a3);

        if(p3 > 0)
        {
          w = (tmax-e0)/tmax ;
          if(p3 > nmaxCont2)
          {
            dp3 = G4float(p3) ;
            corrfac = dp3/G4float(nmaxCont2) ;
            p3 = nmaxCont2 ;
          }
          else
            corrfac = 1. ;

          for(G4int i=0; i<p3; i++) loss += 1./(1.-w*G4UniformRand()) ;
          loss *= e0*corrfac ;  
        }        
      }
    }
    
  else                              // not so small Step
    {
      // excitation type 1
      if(a1>alim)
      {
        siga=sqrt(a1) ;
        p1 = G4std::max(0,int(G4RandGauss::shoot(a1,siga)+0.5));
      }
      else
       p1 = G4Poisson(a1);

      // excitation type 2
      if(a2>alim)
      {
        siga=sqrt(a2) ;
        p2 = G4std::max(0,int(G4RandGauss::shoot(a2,siga)+0.5));
      }
      else
        p2 = G4Poisson(a2);

      loss = p1*e1Fluct+p2*e2Fluct;
 
      // smearing to avoid unphysical peaks
      if(p2 > 0)
        loss += (1.-2.*G4UniformRand())*e2Fluct;   
      else if (loss>0.)
        loss += (1.-2.*G4UniformRand())*e1Fluct;   

      // ionisation .......................................
     if(a3 > 0.)
     {
      if(a3>alim)
      {
        siga=sqrt(a3) ;
        p3 = G4std::max(0,int(G4RandGauss::shoot(a3,siga)+0.5));
      }
      else
        p3 = G4Poisson(a3);

      lossc = 0.;
      if(p3 > 0)
      {
        na = 0.; 
        alfa = 1.;
        if (p3 > nmaxCont2)
        {
          dp3        = G4float(p3);
          rfac       = dp3/(G4float(nmaxCont2)+dp3);
          namean     = G4float(p3)*rfac;
          sa         = G4float(nmaxCont1)*rfac;
          na         = G4RandGauss::shoot(namean,sa);
          if (na > 0.)
          {
            alfa   = w1*G4float(nmaxCont2+p3)/
                    (w1*G4float(nmaxCont2)+G4float(p3));
            alfa1  = alfa*log(alfa)/(alfa-1.);
            ea     = na*ipotFluct*alfa1;
            sea    = ipotFluct*sqrt(na*(alfa-alfa1*alfa1));
            lossc += G4RandGauss::shoot(ea,sea);
          }
        }

        nb = G4int(G4float(p3)-na);
        if (nb > 0)
        {
          w2 = alfa*ipotFluct;
          w  = (tmax-w2)/tmax;      
          for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
        }
      }        
      loss += lossc;  
     }
    } 

  meanLoss = loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4UniversalFluctuation
//
// Author:        Laszlo Urban
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
// 03-10-05 energy dependent rate -> cut dependence of the
//          distribution is much weaker (L.Urban)
// 17-10-05 correction for very small loss (L.Urban)
// 20-03-07 'GLANDZ' part rewritten completely, no 'very small loss'
//          regime any more (L.Urban)
// 03-04-07 correction to get better width of eloss distr.(L.Urban)
// 13-07-07 add protection for very small step or low-density material (VI)
// 19-03-09 new width correction (does not depend on previous steps) (L.Urban)
// 20-03-09 modification in the width correction (L.Urban)
// 14-06-10 fixed tail distribution - do not use uniform function (L.Urban)
// 08-08-10 width correction algorithm has bee modified -->
//          better results for thin targets (L.Urban)
// 06-02-11 correction for very small losses (L.Urban)
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UniversalFluctuation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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
  nmaxCont(16.),
  rate(0.55),
  fw(4.)
{
  lastMaterial = 0;

  particleMass = chargeSquare = ipotFluct = electronDensity = f1Fluct = f2Fluct 
    = e1Fluct = e2Fluct = e1LogFluct = e2LogFluct = ipotLogFluct = e0 = esmall 
    = e1 = e2 = 0;

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
						    G4double& averageLoss)
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

  if(!particle) { InitialiseMe(dp->GetDefinition()); }
  
  G4double tau   = tkin/particleMass;
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
      siga = sqrt((1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
		  * electronDensity * chargeSquare);

     
      G4double sn = meanLoss/siga;
  
      // thick target case 
      if (sn >= 2.0) {

	G4double twomeanLoss = meanLoss + meanLoss;
	do {
	  loss = G4RandGauss::shoot(meanLoss,siga);
	} while (0.0 > loss || twomeanLoss < loss);

	// Gamma distribution
      } else {

	G4double neff = sn*sn;
	loss = meanLoss*CLHEP::RandGamma::shoot(neff,1.0)/neff;
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

  G4double losstot = 0.;
  G4int    nstep = 1;
  if(meanLoss < 25.*ipotFluct)
    {
      if(G4UniformRand() < 0.04*meanLoss/ipotFluct)
	{ nstep = 1; }
      else
	{ 
	  nstep = 2;
	  meanLoss *= 0.5; 
	}
    }

  for (G4int istep=0; istep < nstep; ++istep) {
    
    loss = 0.;

    G4double a1 = 0. , a2 = 0., a3 = 0. ;

    if(tmax > ipotFluct) {
      G4double w2 = log(2.*electron_mass_c2*beta2*gam2)-beta2;

      if(w2 > ipotLogFluct)  {
	G4double C = meanLoss*(1.-rate)/(w2-ipotLogFluct);
	a1 = C*f1Fluct*(w2-e1LogFluct)/e1Fluct;
	if(w2 > e2LogFluct) {
	  a2 = C*f2Fluct*(w2-e2LogFluct)/e2Fluct;
	}
	if(a1 < nmaxCont) { 
	  //small energy loss
	  G4double sa1 = sqrt(a1);
	  if(G4UniformRand() < exp(-sa1))
	    {
	      e1 = esmall;
	      a1 = meanLoss*(1.-rate)/e1;
	      a2 = 0.;
	      e2 = e2Fluct;
	    } 
	  else
	    {
	      a1 = sa1 ;    
	      e1 = sa1*e1Fluct;
	      e2 = e2Fluct;
	    }

	} else {
	    //not small energy loss
	    //correction to get better fwhm value
	    a1 /= fw;
	    e1 = fw*e1Fluct;
	    e2 = e2Fluct;
	}
      }   
    }

    G4double w1 = tmax/e0;
    if(tmax > e0) {
      a3 = rate*meanLoss*(tmax-e0)/(e0*tmax*log(w1));
    }
    //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont  
    G4double emean = 0.;
    G4double sig2e = 0., sige = 0.;
    G4double p1 = 0., p2 = 0., p3 = 0.;
 
    // excitation of type 1
    if(a1 > nmaxCont)
      {
	emean += a1*e1;
	sig2e += a1*e1*e1;
      }
    else if(a1 > 0.)
      {
	p1 = G4double(G4Poisson(a1));
	loss += p1*e1;
	if(p1 > 0.) {
	  loss += (1.-2.*G4UniformRand())*e1;
	}
      }


    // excitation of type 2
    if(a2 > nmaxCont)
      {
	emean += a2*e2;
	sig2e += a2*e2*e2;
      }
    else if(a2 > 0.)
      {
	p2 = G4double(G4Poisson(a2));
	loss += p2*e2;
	if(p2 > 0.) 
	  loss += (1.-2.*G4UniformRand())*e2;
      }

    if(emean > 0.)
      {
	sige   = sqrt(sig2e);
	loss += max(0.,G4RandGauss::shoot(emean,sige));
      }

    // ionisation 
    G4double lossc = 0.;
    if(a3 > 0.) {
      emean = 0.;
      sig2e = 0.;
      sige = 0.;
      p3 = a3;
      G4double alfa = 1.;
      if(a3 > nmaxCont)
	{
	  alfa            = w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
	  G4double alfa1  = alfa*log(alfa)/(alfa-1.);
	  G4double namean = a3*w1*(alfa-1.)/((w1-1.)*alfa);
	  emean          += namean*e0*alfa1;
	  sig2e          += e0*e0*namean*(alfa-alfa1*alfa1);
	  p3              = a3-namean;
	}

      G4double w2 = alfa*e0;
      G4double w  = (tmax-w2)/tmax;
      G4int nb = G4Poisson(p3);
      if(nb > 0) {
	for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
      }
    }

    if(emean > 0.)
      {
	sige   = sqrt(sig2e);
	lossc += max(0.,G4RandGauss::shoot(emean,sige));
      }

    loss += lossc;

    losstot += loss;
  }
  //G4cout << "Vavilov: " << losstot << "  Nstep= " << nstep << G4endl;
              
  return losstot;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4UniversalFluctuation::Dispersion(
                          const G4Material* material,
                          const G4DynamicParticle* dp,
 				G4double& tmax,
			        G4double& length)
{
  if(!particle) { InitialiseMe(dp->GetDefinition()); }

  electronDensity = material->GetElectronDensity();

  G4double gam   = (dp->GetKineticEnergy())/particleMass + 1.0;
  G4double beta2 = 1.0 - 1.0/(gam*gam);

  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                 * electronDensity * chargeSquare;

  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4UniversalFluctuation::SetParticleAndCharge(const G4ParticleDefinition* part,
					     G4double q2)
{
  if(part != particle) {
    particle       = part;
    particleMass   = part->GetPDGMass();
  }
  chargeSquare = q2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

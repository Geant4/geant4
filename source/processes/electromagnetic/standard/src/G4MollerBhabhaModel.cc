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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4MollerBhabhaModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
// Class Description: 
//
// Implementation of energy loss and delta-electron production by e+/e-
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MollerBhabhaModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MollerBhabhaModel::G4MollerBhabhaModel(const G4ParticleDefinition* p) 
  : G4VEmModel(),
  particle(0),
  highKinEnergy(100.*TeV),
  lowKinEnergy(0.1*keV),
  twoln10(2.0*log(10.0)),
  lowLimit(0.2*keV),
  isElectron(true)
{
  if(p) SetParticle(p);
  theElectron = G4Electron::Electron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MollerBhabhaModel::~G4MollerBhabhaModel() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MollerBhabhaModel::SetParticle(const G4ParticleDefinition* p) 
{
  particle = p;
  if(p != theElectron) isElectron = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MollerBhabhaModel::HighEnergyLimit(const G4ParticleDefinition* p,
                                            const G4Material*) 
{
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4double G4MollerBhabhaModel::LowEnergyLimit(const G4ParticleDefinition* p,
                                           const G4Material*) 
{
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MollerBhabhaModel::MinEnergyCut(const G4ParticleDefinition* p,
                                         const G4Material* material) 
{
  return material->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4bool G4MollerBhabhaModel::IsInCharge(const G4ParticleDefinition* p,
			             const G4Material*) 
{
  return (p == theElectron || p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MollerBhabhaModel::ComputeDEDX(const G4Material* material,
                                        const G4ParticleDefinition* p,
                                              G4double kineticEnergy,
                                              G4double cutEnergy) 
{
  if(!particle) SetParticle(p);
  // calculate the dE/dx due to the ionization by Seltzer-Berger formula

  G4double electronDensity = material->GetElectronDensity();
  G4double Zeff  = electronDensity/material->GetTotNbOfAtomsPerVolume();
  G4double th    = 0.25*sqrt(Zeff)*keV;
  G4double tkin  = kineticEnergy;
  if (kineticEnergy < th) tkin = th;
 
  G4double tau   = tkin/electron_mass_c2;
  G4double gam   = tau + 1.0;    
  G4double gamma2= gam*gam;    
  G4double beta2 = 1. - 1./gamma2;
  G4double bg2   = beta2*gamma2;

  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  eexc          /= electron_mass_c2;
  G4double eexc2 = eexc*eexc; 

  G4double dedx;
  G4double d = G4std::min(cutEnergy, MaxSecondaryEnergy(p, tkin))/electron_mass_c2;

  // electron
  if (isElectron) {

    dedx = log(2.0*(tau + 2.0)/eexc2) - 1.0 - beta2
         + log((tau-d)*d) + tau/(tau-d)
         + (0.5*d*d + (2.0*tau + 1.)*log(1. - d/tau))/gamma2;
   
  //positron
  } else {

    G4double d2 = d*d*0.5;
    G4double d3 = d2*d/1.5;
    G4double d4 = d3*d*3.75;
    G4double y  = 1.0/(1.0 + gam);
    dedx = log(2.0*(tau + 2.0)/eexc2) + log(tau*d)
         - beta2*(tau + 2.0*d - y*(3.0*d2 
         + y*(d - d3 + y*(d2 - tau*d3 + d4))))/tau;
  } 

  //density correction 
  G4double cden  = material->GetIonisation()->GetCdensity();
  G4double mden  = material->GetIonisation()->GetMdensity();
  G4double aden  = material->GetIonisation()->GetAdensity();
  G4double x0den = material->GetIonisation()->GetX0density();
  G4double x1den = material->GetIonisation()->GetX1density();
  G4double x     = log(bg2)/twoln10;

  if (x >= x0den) {
    dedx -= twoln10*x - cden;
    if (x < x1den) dedx -= aden*pow(x1den-x, mden);
  } 

  // now you can compute the total ionization loss
  dedx *= twopi_mc2_rcl2*electronDensity/beta2;
  if (dedx < 0.0) dedx = 0.0;

  // lowenergy extrapolation

  if (kineticEnergy < tkin) {

    if (kineticEnergy >= lowLimit) dedx *= sqrt(kineticEnergy/tkin);
    else                           dedx *= sqrt(kineticEnergy*tkin)/lowLimit;

  }
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MollerBhabhaModel::CrossSection(const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy) 
{
  if(!particle) SetParticle(p);
  G4double cross = 0.0;
  G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  tmax = G4std::min(maxEnergy, tmax);

  if(cutEnergy < tmax) {
    
    G4double xmin  = cutEnergy/kineticEnergy;
    G4double xmax  = tmax/kineticEnergy;
    G4double gam   = kineticEnergy/electron_mass_c2 + 1.0;
    G4double gamma2= gam*gam;
    G4double beta2 = 1.0 - 1.0/gamma2;

    //Moller (e-e-) scattering
    if (isElectron) {

      G4double g     = (2.0*gam - 1.0)/gamma2;
      cross = ((xmax - xmin)*(1.0 - g + 1.0/(xmin*xmax) 
			      + 1.0/((1.0-xmin)*(1.0 - xmax))) 
            - g*log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;

    //Bhabha (e+e-) scattering
    } else {

      G4double y   = 1.0/(1.0 + gam);
      G4double y2  = y*y;
      G4double y12 = 1.0 - 2.0*y;
      G4double b1  = 2.0 - y2;
      G4double b2  = y12*(3.0 + y2);
      G4double y122= y12*y12;
      G4double b4  = y122*y12; 
      G4double b3  = b4 + y122;

      cross = (xmax - xmin)*(1.0/(beta2*xmin*xmax) + b2 
            - 0.5*b3*(xmin + xmax) 
	    + b4*(xmin*xmin + xmin*xmax + xmax*xmax)/3.0)
            - b1*log(xmax/xmin);
    }
    
    cross *= twopi_mc2_rcl2*(material->GetElectronDensity())/kineticEnergy;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4DynamicParticle*>* G4MollerBhabhaModel::SampleSecondary(
                             const G4Material* material,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy) 
{
  G4double tmax = G4std::min(maxEnergy, MaxSecondaryEnergy(dp));
  if(tmin >= tmax) return 0;

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double energy = kineticEnergy + electron_mass_c2;
  G4double totalMomentum = sqrt(kineticEnergy*(energy + electron_mass_c2)); 
  G4double xmin   = tmin/kineticEnergy;
  G4double xmax   = tmax/kineticEnergy;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  G4double beta2  = 1.0 - 1.0/gamma2;
  G4double x, z, q, grej;

  G4ThreeVector momentum = dp->GetMomentum();

  //Moller (e-e-) scattering
  if (isElectron) {

    G4double g = (2.0*gam - 1.0)/gamma2;
    G4double y = 1.0 - xmax;
    grej = 1.0 - g*xmax + xmax*xmax*(1.0 - g + (1.0 - g*y)/(y*y));
    
    do {
      q = G4UniformRand();
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = 1.0 - x;
      z = 1.0 - g*x + x*x*(1.0 - g + (1.0 - g*y)/(y*y));
      if(z > grej) {
        G4cout << "G4MollerBhabhaModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << z << " for x= " << x
               << " e-e- scattering"
               << G4endl; 
      }
    } while(grej * G4UniformRand() > z);

  //Bhabha (e+e-) scattering
  } else {

    G4double y   = 1.0/(1.0 + gam);
    G4double y2  = y*y;
    G4double y12 = 1.0 - 2.0*y;
    G4double b1  = 2.0 - y2;
    G4double b2  = y12*(3.0 + y2);
    G4double y122= y12*y12;
    G4double b4  = y122*y12; 
    G4double b3  = b4 + y122;

    y     = xmax*xmax;
    grej  = -xmin*b1;
    grej += y*b2;
    grej -= xmin*xmin*xmin*b3;
    grej += y*y*b4;
    grej *= beta2;
    grej += 1.0;
    do {
      q  = G4UniformRand();
      x  = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      z  = -x*b1;
      y  = x*x;
      z += y*b2;
      y *= x;
      z -= y*b3;
      y *= x;
      z += y*b4;
      z *= beta2;
      z += 1.0;
      if(z > grej) {
        G4cout << "G4MollerBhabhaModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << z << " for x= " << x
               << " e+e- scattering"
               << G4endl; 
      }
    } while(grej * G4UniformRand() > z);
  }
  
  G4double deltaKinEnergy = x * kineticEnergy;
    
  G4double deltaMomentum = 
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
                                   (deltaMomentum * totalMomentum);
  G4double sint = sqrt(1.0 - cost*cost);
 
  G4double phi = twopi * G4UniformRand() ; 

  G4ThreeVector deltaDirection(sint*cos(phi),sint*sin(phi), cost) ;
  deltaDirection.rotateUz(momentum.unit());

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = new G4DynamicParticle(theElectron,
                                                   deltaKinEnergy,
                                                   deltaDirection);
  G4std::vector<G4DynamicParticle*>* vdp = new G4std::vector<G4DynamicParticle*>;
  vdp->push_back(delta);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

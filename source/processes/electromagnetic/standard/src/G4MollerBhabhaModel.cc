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
// File name:   G4MollerBhabhaModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 13-11-02 Minor fix - use normalised direction (V.Ivanchenko)
// 04-12-02 Change G4DynamicParticle constructor in PostStepDoIt (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 25-07-05 Add protection in calculation of recoil direction for the case 
//          of complete energy transfer from e+ to e- (V.Ivanchenko)
// 06-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 15-05-06 Fix MinEnergyCut (V.Ivanchenko)
//
//
// Class Description:
//
// Implementation of energy loss and delta-electron production by e+/e-
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MollerBhabhaModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4Log.hh"
#include "G4DeltaAngle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MollerBhabhaModel::G4MollerBhabhaModel(const G4ParticleDefinition* p,
                                         const G4String& nam)
  : G4VEmModel(nam),
    particle(nullptr),
    isElectron(true),
    twoln10(2.0*G4Log(10.0)),
    lowLimit(0.02*keV),
    isInitialised(false)
{
  theElectron = G4Electron::Electron();
  if(nullptr != p) { SetParticle(p); }
  fParticleChange = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MollerBhabhaModel::~G4MollerBhabhaModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MollerBhabhaModel::MaxSecondaryEnergy(const G4ParticleDefinition*,
                                                 G4double kinEnergy) 
{
  G4double tmax = kinEnergy;
  if(isElectron) { tmax *= 0.5; }
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MollerBhabhaModel::Initialise(const G4ParticleDefinition* p,
                                     const G4DataVector&)
{
  if(p != particle) { SetParticle(p); }

  if(isInitialised) { return; }

  isInitialised = true;
  fParticleChange = GetParticleChangeForLoss();
  if(UseAngularGeneratorFlag() && !GetAngularDistribution()) {
    SetAngularDistribution(new G4DeltaAngle());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MollerBhabhaModel::ComputeCrossSectionPerElectron(
         const G4ParticleDefinition* p, G4double kineticEnergy,
	 G4double cutEnergy, G4double maxEnergy)
{
  if(p != particle) { SetParticle(p); }

  G4double cross = 0.0;
  G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  tmax = std::min(maxEnergy, tmax);
  //G4cout << "E= " << kineticEnergy << " cut= " << cutEnergy 
  //         << " Emax= " << tmax << G4endl;
  if(cutEnergy < tmax) {

    G4double xmin  = cutEnergy/kineticEnergy;
    G4double xmax  = tmax/kineticEnergy;
    G4double tau   = kineticEnergy/electron_mass_c2;
    G4double gam   = tau + 1.0;
    G4double gamma2= gam*gam;
    G4double beta2 = tau*(tau + 2)/gamma2;

    //Moller (e-e-) scattering
    if (isElectron) {

      G4double gg = (2.0*gam - 1.0)/gamma2;
      cross = ((xmax - xmin)*(1.0 - gg + 1.0/(xmin*xmax)
                              + 1.0/((1.0-xmin)*(1.0 - xmax)))
            - gg*G4Log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;

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
            - b1*G4Log(xmax/xmin);
    }

    cross *= twopi_mc2_rcl2/kineticEnergy;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MollerBhabhaModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  return Z*ComputeCrossSectionPerElectron(p,kineticEnergy,cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MollerBhabhaModel::CrossSectionPerVolume(
                                           const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kinEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double eDensity = material->GetElectronDensity();
  return eDensity*ComputeCrossSectionPerElectron(p,kinEnergy,cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MollerBhabhaModel::ComputeDEDXPerVolume(
                                          const G4Material* material,
                                          const G4ParticleDefinition* p,
                                                G4double kineticEnergy,
                                                G4double cut)
{
  if(p != particle) { SetParticle(p); }
  // calculate the dE/dx due to the ionization by Seltzer-Berger formula
  // checl low-energy limit
  G4double electronDensity = material->GetElectronDensity();
  
  G4double Zeff  = material->GetIonisation()->GetZeffective();
  G4double th    = 0.25*sqrt(Zeff)*keV;
  G4double tkin = std::max(kineticEnergy, th);
 
  G4double tau   = tkin/electron_mass_c2;
  G4double gam   = tau + 1.0;
  G4double gamma2= gam*gam;
  G4double bg2   = tau*(tau + 2);
  G4double beta2 = bg2/gamma2;

  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  eexc          /= electron_mass_c2;
  G4double eexc2 = eexc*eexc; 

  G4double d = std::min(cut, MaxSecondaryEnergy(p, tkin))/electron_mass_c2;
  G4double dedx;

  // electron
  if (isElectron) {

    dedx = G4Log(2.0*(tau + 2.0)/eexc2) - 1.0 - beta2
         + G4Log((tau-d)*d) + tau/(tau-d)
         + (0.5*d*d + (2.0*tau + 1.)*G4Log(1. - d/tau))/gamma2;
   
  //positron
  } else {

    G4double d2 = d*d*0.5;
    G4double d3 = d2*d/1.5;
    G4double d4 = d3*d*0.75;
    G4double y  = 1.0/(1.0 + gam);
    dedx = G4Log(2.0*(tau + 2.0)/eexc2) + G4Log(tau*d)
         - beta2*(tau + 2.0*d - y*(3.0*d2 
         + y*(d - d3 + y*(d2 - tau*d3 + d4))))/tau;
  } 

  //density correction 
  G4double x = G4Log(bg2)/twoln10;
  dedx -= material->GetIonisation()->DensityCorrection(x); 

  // now you can compute the total ionization loss
  dedx *= twopi_mc2_rcl2*electronDensity/beta2;
  if (dedx < 0.0) { dedx = 0.0; }
 
  // lowenergy extrapolation
 
  if (kineticEnergy < th) {
    x = kineticEnergy/th;
    if(x > 0.25) { dedx /= sqrt(x); }
    else { dedx *= 1.4*sqrt(x)/(0.1 + x); }
  }
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4MollerBhabhaModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
                                       const G4MaterialCutsCouple* couple,
                                       const G4DynamicParticle* dp,
                                       G4double cutEnergy,
                                       G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  //G4cout << "G4MollerBhabhaModel::SampleSecondaries: E= " << kineticEnergy
  //       << " in " << couple->GetMaterial()->GetName() << G4endl;
  G4double tmax;
  G4double tmin = cutEnergy;  
  if(isElectron) { 
    tmax = 0.5*kineticEnergy; 
  } else {
    tmax = kineticEnergy; 
  }
  if(maxEnergy < tmax) { tmax = maxEnergy; }
  if(tmin >= tmax) { return; }

  G4double energy = kineticEnergy + electron_mass_c2;
  G4double xmin   = tmin/kineticEnergy;
  G4double xmax   = tmax/kineticEnergy;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  G4double beta2  = 1.0 - 1.0/gamma2;
  G4double x, z, grej;
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  G4double rndm[2];

  //Moller (e-e-) scattering
  if (isElectron) {

    G4double gg = (2.0*gam - 1.0)/gamma2;
    G4double y = 1.0 - xmax;
    grej = 1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));

    do {
      rndmEngine->flatArray(2, rndm);
      x = xmin*xmax/(xmin*(1.0 - rndm[0]) + xmax*rndm[0]);
      y = 1.0 - x;
      z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
      /*
      if(z > grej) {
        G4cout << "G4MollerBhabhaModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << z << " for x= " << x
               << " e-e- scattering"
               << G4endl;
      }
      */
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(grej * rndm[1] > z);

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

    y    = xmax*xmax;
    grej = 1.0 + (y*y*b4 - xmin*xmin*xmin*b3 + y*b2 - xmin*b1)*beta2; 
    do {
      rndmEngine->flatArray(2, rndm);
      x = xmin*xmax/(xmin*(1.0 - rndm[0]) + xmax*rndm[0]);
      y = x*x;
      z = 1.0 + (y*y*b4 - x*y*b3 + y*b2 - x*b1)*beta2; 
      /*
      if(z > grej) {
        G4cout << "G4MollerBhabhaModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << z << " for x= " << x
               << " e+e- scattering"
               << G4endl;
      }
      */
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(grej * rndm[1] > z);
  }

  G4double deltaKinEnergy = x * kineticEnergy;

  G4ThreeVector deltaDirection;

  if(UseAngularGeneratorFlag()) {
    const G4Material* mat =  couple->GetMaterial();
    G4int Z = SelectRandomAtomNumber(mat);

    deltaDirection = 
      GetAngularDistribution()->SampleDirection(dp, deltaKinEnergy, Z, mat);

  } else {
 
    G4double deltaMomentum =
      sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
    G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
      (deltaMomentum * dp->GetTotalMomentum());
    if(cost > 1.0) { cost = 1.0; }
    G4double sint = sqrt((1.0 - cost)*(1.0 + cost));

    G4double phi = twopi * rndmEngine->flat() ;

    deltaDirection.set(sint*cos(phi),sint*sin(phi), cost) ;
    deltaDirection.rotateUz(dp->GetMomentumDirection());
  }  

  // create G4DynamicParticle object for delta ray
  auto delta = new G4DynamicParticle(theElectron,deltaDirection,deltaKinEnergy);
  vdp->push_back(delta);

  // primary change
  kineticEnergy -= deltaKinEnergy;
  G4ThreeVector finalP = dp->GetMomentum() - delta->GetMomentum();
  finalP               = finalP.unit();

  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(finalP);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

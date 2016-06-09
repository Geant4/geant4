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
// $Id: G4LowEnergyBraggModel.cc,v 1.6 2006/06/29 19:40:11 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4LowEnergyBraggModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 17.04.2004
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss and delta-electron production by
// slow charged heavy particles

// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LowEnergyBraggModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4hParametrisedLossModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LowEnergyBraggModel::G4LowEnergyBraggModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),
  particle(0),
  parameterization(0),
  isIon(false)
{
  if(p) SetParticle(p);
  highKinEnergy = 10.0*MeV;
  lowKinEnergy = 0.0*MeV;
  lowestKinEnergy = 1.0*keV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LowEnergyBraggModel::~G4LowEnergyBraggModel()
{
  delete parameterization;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBraggModel::SetParticle(const G4ParticleDefinition* p)
{
  if(particle != p) {
    particle = p;
    mass = particle->GetPDGMass();
    spin = particle->GetPDGSpin();
    G4double q = particle->GetPDGCharge()/eplus;
    chargeSquare = q*q;
    massRate     = mass/proton_mass_c2;
//    highKinEnergy *= massRate;
    ratio = electron_mass_c2/mass;
    if(particle->GetParticleName() == "GenericIon") isIon = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyBraggModel::HighEnergyLimit(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyBraggModel::LowEnergyLimit(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBraggModel::SetParameterization(G4hParametrisedLossModel* p)
{
  delete parameterization;
  parameterization = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyBraggModel::MinEnergyCut(const G4ParticleDefinition*,
                                    const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4LowEnergyBraggModel::IsInCharge(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return (p->GetPDGCharge() != 0.0 && p->GetPDGMass() > 10.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBraggModel::Initialise(const G4ParticleDefinition* p,
                                       const G4DataVector&)
{
  if(!particle) SetParticle(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyBraggModel::ComputeDEDX(const G4MaterialCutsCouple* couple,
                                            const G4ParticleDefinition* p,
                                                  G4double kineticEnergy,
                                                  G4double cutEnergy)
{
  const G4Material* material = couple->GetMaterial();
  G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tkin  = kineticEnergy*massRate;
  G4double dedx;
  if(tkin >= lowestKinEnergy) 
    dedx  = parameterization->TheValue(p, material, kineticEnergy);
  else {
    G4double e = lowestKinEnergy/massRate;
    dedx  = (parameterization->TheValue(p, material, e))*std::sqrt(kineticEnergy/e);
  }

  if (cutEnergy < tmax) {

    G4double tau   = kineticEnergy/mass;
    G4double gam   = tau + 1.0;
    G4double bg2   = tau * (tau+2.0);
    G4double beta2 = bg2/(gam*gam);
    G4double x     = cutEnergy/tmax;

    G4double delta = std::log(x) - (1.0 - x)*beta2;
    if(spin == 0.5) {
      G4double y = tmax/(kineticEnergy + mass);
      delta += 0.25*y*y*(1.0 - x*x);
    }

    delta *= twopi_mc2_rcl2*(material->GetElectronDensity())/beta2;
    dedx  -= delta;
  }

  // now compute the total ionization loss

  if (dedx < 0.0) dedx = 0.0 ;

  dedx *= chargeSquare;

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyBraggModel::CrossSection(const G4MaterialCutsCouple* couple,
                                    const G4ParticleDefinition* p,
                                          G4double kineticEnergy,
                                          G4double cutEnergy,
                                          G4double maxKinEnergy)
{

  G4double cross     = 0.0;
  G4double tmax      = MaxSecondaryEnergy(p, kineticEnergy);
  G4double maxEnergy = std::min(tmax,maxKinEnergy);
  if(cutEnergy < tmax) {

    G4double energy  = kineticEnergy + mass;
    G4double energy2 = energy*energy;
    G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;
    cross = 1.0/cutEnergy - 1.0/maxEnergy - beta2*std::log(maxEnergy/cutEnergy)/tmax;

// +term for spin=1/2 particle
    if( 0.5 == spin ) {
      cross        +=  0.5 * (tmax - cutEnergy) / energy2;
    }
    cross *= twopi_mc2_rcl2*chargeSquare*
             (couple->GetMaterial()->GetElectronDensity())/beta2;
  }
 //   G4cout << "BR: e= " << kineticEnergy << " tmin= " << cutEnergy << " tmax= " << tmax
 //        << " cross= " << cross << G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* G4LowEnergyBraggModel::SampleSecondary(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{
  G4double tmax = MaxSecondaryEnergy(dp);
  G4double xmax = std::min(tmax, maxEnergy);
  G4double xmin = std::min(xmax,tmin);

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double energy  = kineticEnergy + mass;
  G4double energy2 = energy*energy;
  G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;
  G4double grej    = 1.0;
  G4double deltaKinEnergy, f;

  G4ThreeVector momentum = dp->GetMomentumDirection();

  // sampling follows ...
  do {
    G4double q = G4UniformRand();
    deltaKinEnergy = xmin*xmax/(xmin*(1.0 - q) + xmax*q);

    f = 1.0 - beta2*deltaKinEnergy/tmax;

    if(f > grej) {
        G4cout << "G4LowEnergyBraggModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for e= " << deltaKinEnergy
               << G4endl;
    }

  } while( grej*G4UniformRand() >= f );

  G4double deltaMomentum =
           std::sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double totMomentum = std::sqrt(energy2 - mass*mass);
  G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
                                   (deltaMomentum * totMomentum);
  G4double sint = std::sqrt(1.0 - cost*cost);

  G4double phi = twopi * G4UniformRand() ;

  G4ThreeVector deltaDirection(sint*std::cos(phi),sint*std::sin(phi), cost) ;
  deltaDirection.rotateUz(momentum);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = new G4DynamicParticle(G4Electron::Electron(),
                                                   deltaDirection,deltaKinEnergy);

  return delta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>* G4LowEnergyBraggModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{
  std::vector<G4DynamicParticle*>* vdp = new std::vector<G4DynamicParticle*>;
  G4DynamicParticle* delta = SampleSecondary(couple, dp, tmin, maxEnergy);
  vdp->push_back(delta);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


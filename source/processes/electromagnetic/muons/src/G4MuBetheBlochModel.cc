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
// GEANT4 Class header file
//
//
// File name:     G4MuBetheBlochModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 09.08.2002
//
// Modifications:
//
// 04-12-02 Fix problem of G4DynamicParticle constructor (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 10-02-04 Calculation of radiative corrections using R.Kokoulin model (V.I)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 12-04-05 Add usage of G4EmCorrections (V.Ivanchenko)
// 13-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
//

//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuBetheBlochModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4EmCorrections.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

G4double G4MuBetheBlochModel::xgi[]={ 0.0199, 0.1017, 0.2372, 0.4083, 0.5917,
                                      0.7628, 0.8983, 0.9801 };
                                      
G4double G4MuBetheBlochModel::wgi[]={ 0.0506, 0.1112, 0.1569, 0.1813, 0.1813,
                                      0.1569, 0.1112, 0.0506 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuBetheBlochModel::G4MuBetheBlochModel(const G4ParticleDefinition* p,
                                         const G4String& nam)
  : G4VEmModel(nam),
    limitKinEnergy(100.*CLHEP::keV),
    logLimitKinEnergy(G4Log(limitKinEnergy)),
    twoln10(2.0*G4Log(10.0)),
    alphaprime(CLHEP::fine_structure_const/CLHEP::twopi)
{
  theElectron = G4Electron::Electron();
  corr = G4LossTableManager::Instance()->EmCorrections();
  if(nullptr != p) { SetParticle(p); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBetheBlochModel::MinEnergyCut(const G4ParticleDefinition*,
                                           const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4MuBetheBlochModel::MaxSecondaryEnergy(const G4ParticleDefinition*,
                                                 G4double kinEnergy) 
{
  G4double tau  = kinEnergy/mass;
  G4double tmax = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.) /
                  (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuBetheBlochModel::SetParticle(const G4ParticleDefinition* p)
{
  if(nullptr == particle) {
    particle = p;
    mass = particle->GetPDGMass();
    massSquare = mass*mass;
    ratio = CLHEP::electron_mass_c2/mass;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuBetheBlochModel::Initialise(const G4ParticleDefinition* p,
                                     const G4DataVector&)
{
  SetParticle(p);
  if(nullptr == fParticleChange) { 
    fParticleChange = GetParticleChangeForLoss();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBetheBlochModel::ComputeCrossSectionPerElectron(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxKinEnergy)
{
  G4double cross = 0.0;
  G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  G4double maxEnergy = std::min(tmax, maxKinEnergy);
  if(cutEnergy < maxEnergy) {

    G4double totEnergy = kineticEnergy + mass;
    G4double energy2 = totEnergy*totEnergy;
    G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;

    cross = 1.0/cutEnergy - 1.0/maxEnergy - 
      beta2*G4Log(maxEnergy/cutEnergy)/tmax +
      0.5*(maxEnergy - cutEnergy)/energy2;

    // radiative corrections of R. Kokoulin
    if (maxEnergy > limitKinEnergy) {

      G4double logtmax = G4Log(maxEnergy);
      G4double logtmin = G4Log(std::max(cutEnergy,limitKinEnergy));
      G4double logstep = logtmax - logtmin;
      G4double dcross  = 0.0;

      for (G4int ll=0; ll<8; ++ll) {
        G4double ep = G4Exp(logtmin + xgi[ll]*logstep);
        G4double a1 = G4Log(1.0 + 2.0*ep/CLHEP::electron_mass_c2);
        G4double a3 = G4Log(4.0*totEnergy*(totEnergy - ep)/massSquare);
        dcross += wgi[ll]*(1.0/ep - beta2/tmax + 0.5*ep/energy2)*a1*(a3 - a1);
      }
      cross += dcross*logstep*alphaprime;
    }
    cross *= CLHEP::twopi_mc2_rcl2/beta2;
  }
  //  G4cout << "tmin= " << cutEnergy << " tmax= " << tmax
  //         << " cross= " << cross << G4endl;  
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBetheBlochModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double cross = Z*ComputeCrossSectionPerElectron
                                         (p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBetheBlochModel::CrossSectionPerVolume(
                                           const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double eDensity = material->GetElectronDensity();
  G4double cross = eDensity*ComputeCrossSectionPerElectron
                                         (p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBetheBlochModel::ComputeDEDXPerVolume(const G4Material* material,
                                                  const G4ParticleDefinition* p,
                                                  G4double kineticEnergy,
                                                  G4double cut)
{
  G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tau   = kineticEnergy/mass;
  G4double cutEnergy = std::min(cut, tmax);
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau+2.0);
  G4double beta2 = bg2/(gam*gam);

  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc2 = eexc*eexc;

  G4double eDensity = material->GetElectronDensity();

  G4double dedx = G4Log(2.0*CLHEP::electron_mass_c2*bg2*cutEnergy/eexc2)
                 -(1.0 + cutEnergy/tmax)*beta2;

  G4double totEnergy = kineticEnergy + mass;
  G4double del = 0.5*cutEnergy/totEnergy;
  dedx += del*del;

  // density correction
  G4double x = G4Log(bg2)/twoln10;
  dedx -= material->GetIonisation()->DensityCorrection(x);

  // shell correction
  dedx -= 2.0*corr->ShellCorrection(p,material,kineticEnergy);
  dedx = std::max(dedx, 0.0);

  // radiative corrections of R. Kokoulin
  if (cutEnergy > limitKinEnergy) {

    G4double logtmax = G4Log(cutEnergy);
    G4double logstep = logtmax - logLimitKinEnergy;
    G4double dloss = 0.0;
    G4double ftot2= 0.5/(totEnergy*totEnergy);

    for (G4int ll=0; ll<8; ++ll) {
      G4double ep = G4Exp(logLimitKinEnergy + xgi[ll]*logstep);
      G4double a1 = G4Log(1.0 + 2.0*ep/CLHEP::electron_mass_c2);
      G4double a3 = G4Log(4.0*totEnergy*(totEnergy - ep)/massSquare);
      dloss += wgi[ll]*(1.0 - beta2*ep/tmax + ep*ep*ftot2)*a1*(a3 - a1);
    }
    dedx += dloss*logstep*alphaprime;
  }

  dedx *= CLHEP::twopi_mc2_rcl2*eDensity/beta2;

  //High order corrections
  dedx += corr->HighOrderCorrections(p,material,kineticEnergy,cutEnergy);
  dedx = std::max(dedx, 0.);
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuBetheBlochModel::SampleSecondaries(
                          std::vector<G4DynamicParticle*>* vdp,
			  const G4MaterialCutsCouple*,
			  const G4DynamicParticle* dp,
			  G4double minKinEnergy,
			  G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double tmax = MaxSecondaryEnergy(dp->GetDefinition(), kineticEnergy);
  G4double maxKinEnergy = std::min(maxEnergy, tmax);
  if(minKinEnergy >= maxKinEnergy) { return; }

  G4double totEnergy = kineticEnergy + mass;
  G4double etot2 = totEnergy*totEnergy;
  G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*mass)/etot2;
 
  G4double grej  = 1.;
  if(tmax > limitKinEnergy) {
    G4double a0 = G4Log(2.*totEnergy/mass);
    grej += alphaprime*a0*a0;
  }

  G4double tkin, f;

  // sampling follows ...
  do {
    G4double q = G4UniformRand();
    tkin = minKinEnergy*maxKinEnergy/(minKinEnergy*(1.0 - q) + maxKinEnergy*q);
    f = 1.0 - beta2*tkin/tmax + 0.5*tkin*tkin/etot2;

    if(tkin > limitKinEnergy) {
      G4double a1 = G4Log(1.0 + 2.0*tkin/CLHEP::electron_mass_c2);
      G4double a3 = G4Log(4.0*totEnergy*(totEnergy - tkin)/massSquare);
      f *= (1. + alphaprime*a1*(a3 - a1));
    }

    if(f > grej) {
        G4cout << "G4MuBetheBlochModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for edelta= " << tkin
               << " tmin= " << minKinEnergy << " max= " << maxKinEnergy
               << G4endl;
    }
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while( grej*G4UniformRand() > f );

  G4double deltaMomentum =
    std::sqrt(tkin * (tkin + 2.0*CLHEP::electron_mass_c2));
  G4double totalMomentum = totEnergy*std::sqrt(beta2);
  G4double cost = tkin * (totEnergy + CLHEP::electron_mass_c2) /
                                   (deltaMomentum * totalMomentum);

  G4double sint = std::sqrt(1.0 - cost*cost);
  G4double phi = CLHEP::twopi * G4UniformRand();
  G4ThreeVector deltaDirection(sint*std::cos(phi), sint*std::sin(phi), cost);
  G4ThreeVector direction = dp->GetMomentumDirection();
  deltaDirection.rotateUz(direction);

  // primary change
  kineticEnergy -= tkin;
  G4ThreeVector dir = totalMomentum*direction - deltaMomentum*deltaDirection;
  direction = dir.unit();
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(direction);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta =
    new G4DynamicParticle(theElectron, deltaDirection, tkin);
  vdp->push_back(delta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

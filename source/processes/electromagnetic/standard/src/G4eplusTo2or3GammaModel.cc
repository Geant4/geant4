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
// File name:   G4eplusTo2or3GammaModel
//
// Author:      Vladimir Ivanchenko and Omrame Kadri
//
// Creation date: 29.03.2018
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "G4eplusTo2or3GammaModel.hh"
#include "G4eplusTo3GammaOKVIModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4TrackStatus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4DataVector.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4eplusTo2or3GammaModel::fCrossSection   = nullptr;
G4PhysicsVector* G4eplusTo2or3GammaModel::f3GProbability  = nullptr;

G4eplusTo2or3GammaModel::G4eplusTo2or3GammaModel()
  : G4VEmModel("eplusTo2or3gamma"),
    fDeltaMin(0.001),
    fDelta(fDeltaMin),
    fGammaTh(CLHEP::MeV)
{
  theGamma = G4Gamma::Gamma();
  fParticleChange = nullptr;
  f3GModel = new G4eplusTo3GammaOKVIModel();
  SetTripletModel(f3GModel);

  // instantiate vectors once
  if (nullptr == fCrossSection) {
    G4double emin = 10*CLHEP::eV;
    G4double emax = 100*CLHEP::TeV;
    G4int nbins = 20*G4lrint(std::log10(emax/emin));
    fCrossSection = new G4PhysicsLogVector(emin, emax, nbins, true);
    f3GProbability = new G4PhysicsLogVector(emin, emax, nbins, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusTo2or3GammaModel::~G4eplusTo2or3GammaModel()
{
  if (IsMaster()) {
    delete fCrossSection;
    delete f3GProbability;
    fCrossSection = nullptr;
    f3GProbability = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusTo2or3GammaModel::Initialise(const G4ParticleDefinition* p,
					  const G4DataVector& cuts)
{
  // here particle change is set for the triplet model
  if (nullptr == fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
  }
  // initialialise 3-gamma model before new run
  f3GModel->Initialise(p, cuts);
  fGammaTh = G4EmParameters::Instance()->LowestTripletEnergy();

  // initialise vectors
  if (IsMaster()) {
    std::size_t num = fCrossSection->GetVectorLength();
    for (std::size_t i=0; i<num; ++i) {
      G4double e = fCrossSection->Energy(i);
      G4double cs2 = ComputeCrossSectionPerElectron(e);
      G4double cs3 = f3GModel->ComputeCrossSectionPerElectron(e);
      cs2 += cs3;
      fCrossSection->PutValue(i, cs2);
      G4double y = (cs2 > 0.0) ? cs3/cs2 : 0.0;
      f3GProbability->PutValue(i, y);
    }
    fCrossSection->FillSecondDerivatives();
    f3GProbability->FillSecondDerivatives();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4eplusTo2or3GammaModel::ComputeCrossSectionPerElectron(G4double kinEnergy)
{
  // Calculates the cross section per electron of annihilation into two 
  // photons from the Heilter formula with the radiation correction to 3 gamma 
  // annihilation channel. (A.A.) rho is changed

  G4double ekin   = std::max(CLHEP::eV, kinEnergy);   
  G4double tau    = ekin/CLHEP::electron_mass_c2;
  G4double gam    = tau + 1.0;
  G4double gamma2 = gam*gam;
  G4double bg2    = tau * (tau+2.0);
  G4double bg     = std::sqrt(bg2);
  G4double rho = (gamma2+4.*gam+1.)*G4Log(gam+bg)/(gamma2-1.) 
    - (gam+3.)/(std::sqrt(gam*gam - 1.));
  G4double eGammaCMS = CLHEP::electron_mass_c2 * std::sqrt(0.5*(tau + 2.0));
  fDelta = std::max(fDeltaMin, fGammaTh/eGammaCMS);
  f3GModel->SetDelta(fDelta);

  static const G4double pir2 =
    CLHEP::pi*CLHEP::classic_electr_radius*CLHEP::classic_electr_radius;
  G4double cross = (pir2*rho + alpha_rcl2*2.*G4Log(fDelta)*rho*rho)/(gam+1.);

  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eplusTo2or3GammaModel::ComputeCrossSectionPerAtom(
                                    const G4ParticleDefinition*,
                                    G4double kineticEnergy, G4double Z,
				    G4double, G4double, G4double)
{
  // Calculates the cross section per atom of annihilation into two photons
  G4double cross = Z*fCrossSection->Value(kineticEnergy);
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eplusTo2or3GammaModel::CrossSectionPerVolume(
					const G4Material* material,
					const G4ParticleDefinition*,
					      G4double kineticEnergy,
					      G4double, G4double)
{
  // Calculates the cross section per volume of annihilation into two photons
  G4double eDensity = material->GetElectronDensity();
  G4double cross = eDensity*fCrossSection->Value(kineticEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Polarisation of gamma according to M.H.L.Pryce and J.C.Ward, 
// Nature 4065 (1947) 435.

void G4eplusTo2or3GammaModel::SampleSecondaries(
				     std::vector<G4DynamicParticle*>* vdp,
                                     const G4MaterialCutsCouple* couple,
                                     const G4DynamicParticle* dp,
                                     G4double, G4double)
{
  // kill primary positron
  fParticleChange->SetProposedKineticEnergy(0.0);
  fParticleChange->ProposeTrackStatus(fStopAndKill);

  // Case at rest not considered anymore
  G4double posiKinEnergy = dp->GetKineticEnergy();
  G4LorentzVector lv(dp->GetMomentum(),
		     posiKinEnergy + 2*CLHEP::electron_mass_c2);
  G4double eGammaCMS = 0.5 * lv.mag();

  if (G4UniformRand() < f3GProbability->Value(posiKinEnergy)) {
    fDelta = std::max(fDeltaMin, fGammaTh/eGammaCMS);
    f3GModel->SetDelta(fDelta);
    f3GModel->SampleSecondaries(vdp, couple, dp);
    return;
  }

  G4ThreeVector dir1 = G4RandomDirection();
  G4double phi = CLHEP::twopi * G4UniformRand();
  G4double cosphi = std::cos(phi);
  G4double sinphi = std::sin(phi);
  G4ThreeVector pol1(cosphi, sinphi, 0.0);
  pol1.rotateUz(dir1);
  G4LorentzVector lv1(eGammaCMS*dir1, eGammaCMS);
  
  G4ThreeVector pol2(-sinphi, cosphi, 0.0);
  pol2.rotateUz(dir1);

  // transformation to lab system
  lv1.boost(lv.boostVector());
  lv -= lv1;

  //!!! boost of polarisation vector is not yet implemented
  
  // use constructors optimal for massless particle
  auto aGamma1 = new G4DynamicParticle(G4Gamma::Gamma(), lv1.vect());
  aGamma1->SetPolarization(pol1);
  auto aGamma2 = new G4DynamicParticle(G4Gamma::Gamma(), lv.vect());
  aGamma2->SetPolarization(pol2);
 
  vdp->push_back(aGamma1);
  vdp->push_back(aGamma2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

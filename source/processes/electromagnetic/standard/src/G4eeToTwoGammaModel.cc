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
// File name:   G4eeToTwoGammaModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 18-04-05 Compute CrossSectionPerVolume (V.Ivanchenko)
// 06-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 29-06-06 Fix problem for zero energy incident positron (V.Ivanchenko) 
// 20-10-06 Add theGamma as a member (V.Ivanchenko)
// 18-01-20 Introduce thermal model of annihilation at rest (J.Allison)
//
//
// Class Description:
//
// Implementation of e+ annihilation into 2 gamma
//
// The secondaries Gamma energies are sampled using the Heitler cross section.
//
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1: The initial electron is assumed free and at rest if atomic PDF 
//         is not defined
//
// Note 2: The annihilation processes producing one or more than two photons are
//         ignored, as negligible compared to the two photons process.

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToTwoGammaModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackStatus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4EmParameters.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToTwoGammaModel::G4eeToTwoGammaModel(const G4ParticleDefinition*,
                                         const G4String& nam)
  : G4VEmModel(nam),
    pi_rcl2(CLHEP::pi*CLHEP::classic_electr_radius*CLHEP::classic_electr_radius)
{
  theGamma = G4Gamma::Gamma();
  fParticleChange = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToTwoGammaModel::~G4eeToTwoGammaModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToTwoGammaModel::Initialise(const G4ParticleDefinition*,
                                     const G4DataVector&)
{
  if (nullptr != fParticleChange) { return; }
  fParticleChange = GetParticleChangeForGamma();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4eeToTwoGammaModel::ComputeCrossSectionPerElectron(G4double kineticEnergy)
{
  // Calculates the cross section per electron of annihilation into two photons
  // from the Heilter formula.

  G4double ekin  = std::max(CLHEP::eV, kineticEnergy);

  G4double tau   = ekin/CLHEP::electron_mass_c2;
  G4double gam   = tau + 1.0;
  G4double gamma2= gam*gam;
  G4double bg2   = tau * (tau+2.0);
  G4double bg    = std::sqrt(bg2);

  G4double cross = pi_rcl2*((gamma2+4*gam+1.)*G4Log(gam+bg) - (gam+3.)*bg)
                 / (bg2*(gam+1.));
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eeToTwoGammaModel::ComputeCrossSectionPerAtom(
                                    const G4ParticleDefinition*,
                                    G4double kineticEnergy, G4double Z,
				    G4double, G4double, G4double)
{
  // Calculates the cross section per atom of annihilation into two photons
  return Z*ComputeCrossSectionPerElectron(kineticEnergy);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eeToTwoGammaModel::CrossSectionPerVolume(
					const G4Material* material,
					const G4ParticleDefinition*,
					      G4double kineticEnergy,
					      G4double, G4double)
{
  // Calculates the cross section per volume of annihilation into two photons
  return material->GetElectronDensity()*ComputeCrossSectionPerElectron(kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Polarisation of gamma according to M.H.L.Pryce and J.C.Ward,
// Nature 4065 (1947) 435.

void G4eeToTwoGammaModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
					    const G4MaterialCutsCouple*,
					    const G4DynamicParticle* dp,
					    G4double,
					    G4double)
{
  // kill primary positron
  fParticleChange->SetProposedKineticEnergy(0.0);
  fParticleChange->ProposeTrackStatus(fStopAndKill);

  // Case at rest not considered anymore inside this model
  G4LorentzVector lv(dp->GetMomentum(),
		     dp->GetKineticEnergy() + 2*CLHEP::electron_mass_c2);
  G4double eGammaCMS = 0.5 * lv.mag();

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

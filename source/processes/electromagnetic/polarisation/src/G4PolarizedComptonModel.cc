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
// -------------------------------------------------------------------
//
// Geant4 Class file
//
// File name:     G4PolarizedComptonModel
//
// Author:        Andreas Schaelicke

#include "G4PolarizedComptonModel.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PhysicalConstants.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizedComptonXS.hh"
#include "G4StokesVector.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedComptonModel::G4PolarizedComptonModel(const G4ParticleDefinition*,
                                                 const G4String& nam)
  : G4KleinNishinaCompton(nullptr, nam)
  , fVerboseLevel(0)
{
  fCrossSectionCalculator = new G4PolarizedComptonXS();
  fBeamPolarization       = G4StokesVector::ZERO;
  fTargetPolarization     = G4StokesVector::ZERO;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedComptonModel::~G4PolarizedComptonModel()
{
  delete fCrossSectionCalculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedComptonModel::ComputeAsymmetryPerAtom(G4double gammaEnergy,
                                                          G4double /*Z*/)
{
  G4double asymmetry = 0.0;

  G4double k0 = gammaEnergy / electron_mass_c2;
  G4double k1 = 1. + 2. * k0;

  asymmetry = -k0;
  asymmetry *=
    (k0 + 1.) * sqr(k1) * G4Log(k1) - 2. * k0 * (5. * sqr(k0) + 4. * k0 + 1.);
  asymmetry /= ((k0 - 2.) * k0 - 2.) * sqr(k1) * G4Log(k1) +
               2. * k0 * (k0 * (k0 + 1.) * (k0 + 8.) + 2.);

  if(asymmetry > 1.)
  {
    G4ExceptionDescription ed;
    ed << "ERROR in G4PolarizedComptonModel::ComputeAsymmetryPerAtom.\n"
       << " asymmetry = " << asymmetry << "\n";
    G4Exception("G4PolarizedComptonModel::ComputeAsymmetryPerAtom", "pol035",
                JustWarning, ed);
  }

  return asymmetry;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedComptonModel::ComputeCrossSectionPerAtom(
  const G4ParticleDefinition* pd, G4double kinEnergy, G4double Z, G4double A,
  G4double cut, G4double emax)
{
  G4double xs = G4KleinNishinaCompton::ComputeCrossSectionPerAtom(
    pd, kinEnergy, Z, A, cut, emax);
  G4double polzz = fBeamPolarization.p3() * fTargetPolarization.z();
  if(polzz > 0.0)
  {
    G4double asym = ComputeAsymmetryPerAtom(kinEnergy, Z);
    xs *= (1. + polzz * asym);
  }
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedComptonModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* fvect, const G4MaterialCutsCouple*,
  const G4DynamicParticle* aDynamicGamma, G4double, G4double)
{
  // do nothing below the threshold
  if(aDynamicGamma->GetKineticEnergy() <= LowEnergyLimit())
  {
    return;
  }

  const G4Track* aTrack       = fParticleChange->GetCurrentTrack();
  G4VPhysicalVolume* aPVolume = aTrack->GetVolume();
  G4LogicalVolume* aLVolume   = aPVolume->GetLogicalVolume();

  if(fVerboseLevel >= 1)
  {
    G4cout << "G4PolarizedComptonModel::SampleSecondaries in "
           << aLVolume->GetName() << G4endl;
  }
  G4PolarizationManager* polarizationManager =
    G4PolarizationManager::GetInstance();

  // obtain polarization of the beam
  fBeamPolarization = G4StokesVector(aDynamicGamma->GetPolarization());
  fBeamPolarization.SetPhoton();

  // obtain polarization of the media
  G4bool targetIsPolarized = polarizationManager->IsPolarized(aLVolume);
  fTargetPolarization = polarizationManager->GetVolumePolarization(aLVolume);

  // if beam is linear polarized or target is transversely polarized
  // determine the angle to x-axis
  // (assumes same PRF as in the polarization definition)
  G4ThreeVector gamDirection0 = aDynamicGamma->GetMomentumDirection();

  // transfer fTargetPolarization
  // into the gamma frame (problem electron is at rest)
  if(targetIsPolarized)
  {
    fTargetPolarization.rotateUz(gamDirection0);
  }
  // The scattered gamma energy is sampled according to
  // Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used
  // (Nuc Phys 20(1960),15).
  // Note : Effects due to binding of atomic electrons are neglected.

  G4double gamEnergy0 = aDynamicGamma->GetKineticEnergy();
  G4double E0_m       = gamEnergy0 / electron_mass_c2;

  // sample the energy rate of the scattered gamma
  G4double epsilon, sint2;
  G4double onecost = 0.0;
  G4double Phi     = 0.0;
  G4double greject = 1.0;
  G4double cosTeta = 1.0;
  G4double sinTeta = 0.0;

  G4double eps0       = 1. / (1. + 2. * E0_m);
  G4double epsilon0sq = eps0 * eps0;
  G4double alpha1     = -G4Log(eps0);
  G4double alpha2     = alpha1 + 0.5 * (1. - epsilon0sq);

  G4double polarization = fBeamPolarization.p3() * fTargetPolarization.p3();

  CLHEP::HepRandomEngine* rndmEngineMod = G4Random::getTheEngine();
  G4int nloop                           = 0;
  G4bool end                            = false;

  G4double rndm[3];

  do
  {
    do
    {
      ++nloop;
      // false interaction if too many iterations
      if(nloop > fLoopLim)
      {
        PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi,
                     "too many iterations");
        return;
      }

      // 3 random numbers to sample scattering
      rndmEngineMod->flatArray(3, rndm);

      if(alpha1 > alpha2 * rndm[0])
      {
        epsilon = G4Exp(-alpha1 * rndm[1]);
      }
      else
      {
        epsilon = std::sqrt(epsilon0sq + (1. - epsilon0sq) * rndm[1]);
      }

      onecost = (1. - epsilon) / (epsilon * E0_m);
      sint2   = onecost * (2. - onecost);

      G4double gdiced = 2. * (1. / epsilon + epsilon);
      G4double gdist  = 1. / epsilon + epsilon - sint2 -
                       polarization * (1. / epsilon - epsilon) * (1. - onecost);

      greject = gdist / gdiced;

      if(greject > 1.0)
      {
        PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi,
                     "theta majoranta wrong");
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(greject < rndm[2]);

    // assuming phi loop successful
    end = true;

    // scattered gamma angles. ( Z - axis along the parent gamma)
    cosTeta = 1. - onecost;
    sinTeta = std::sqrt(sint2);
    do
    {
      ++nloop;

      // 2 random numbers to sample scattering
      rndmEngineMod->flatArray(2, rndm);

      // false interaction if too many iterations
      Phi = twopi * rndm[0];
      if(nloop > fLoopLim)
      {
        PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi,
                     "too many iterations");
        return;
      }

      G4double gdiced = 1. / epsilon + epsilon - sint2 +
                        std::abs(fBeamPolarization.p3()) *
                          (std::abs((1. / epsilon - epsilon) * cosTeta *
                                    fTargetPolarization.p3()) +
                           (1. - epsilon) * sinTeta *
                             (std::sqrt(sqr(fTargetPolarization.p1()) +
                                        sqr(fTargetPolarization.p2())))) +
                        sint2 * (std::sqrt(sqr(fBeamPolarization.p1()) +
                                           sqr(fBeamPolarization.p2())));

      G4double gdist =
        1. / epsilon + epsilon - sint2 +
        fBeamPolarization.p3() *
          ((1. / epsilon - epsilon) * cosTeta * fTargetPolarization.p3() +
           (1. - epsilon) * sinTeta *
             (std::cos(Phi) * fTargetPolarization.p1() +
              std::sin(Phi) * fTargetPolarization.p2())) -
        sint2 * (std::cos(2. * Phi) * fBeamPolarization.p1() +
                 std::sin(2. * Phi) * fBeamPolarization.p2());
      greject = gdist / gdiced;

      if(greject > 1.0)
      {
        PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi,
                     "phi majoranta wrong");
      }

      if(greject < 1.e-3)
      {
        PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi,
                     "phi loop ineffective");
        // restart theta loop
        end = false;
        break;
      }

      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(greject < rndm[1]);
  } while(!end);
  G4double dirx = sinTeta * std::cos(Phi);
  G4double diry = sinTeta * std::sin(Phi);
  G4double dirz = cosTeta;

  // update G4VParticleChange for the scattered gamma
  G4ThreeVector gamDirection1(dirx, diry, dirz);
  gamDirection1.rotateUz(gamDirection0);
  G4double gamEnergy1 = epsilon * gamEnergy0;

  G4double edep = 0.0;
  if(gamEnergy1 > lowestSecondaryEnergy)
  {
    fParticleChange->ProposeMomentumDirection(gamDirection1);
    fParticleChange->SetProposedKineticEnergy(gamEnergy1);
  }
  else
  {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    edep = gamEnergy1;
  }

  // calculate Stokes vector of final state photon and electron
  G4ThreeVector nInteractionFrame =
    G4PolarizationHelper::GetFrame(gamDirection1, gamDirection0);

  // transfer fBeamPolarization and fTargetPolarization
  // into the interaction frame (note electron is in gamma frame)
  if(fVerboseLevel >= 1)
  {
    G4cout << "========================================" << G4endl;
    G4cout << " nInteractionFrame = " << nInteractionFrame << G4endl;
    G4cout << " GammaDirection0 = " << gamDirection0 << G4endl;
    G4cout << " gammaPolarization = " << fBeamPolarization << G4endl;
    G4cout << " electronPolarization = " << fTargetPolarization << G4endl;
  }

  fBeamPolarization.InvRotateAz(nInteractionFrame, gamDirection0);
  fTargetPolarization.InvRotateAz(nInteractionFrame, gamDirection0);

  if(fVerboseLevel >= 1)
  {
    G4cout << "----------------------------------------" << G4endl;
    G4cout << " gammaPolarization = " << fBeamPolarization << G4endl;
    G4cout << " electronPolarization = " << fTargetPolarization << G4endl;
    G4cout << "----------------------------------------" << G4endl;
  }

  // initialize the polarization transfer matrix
  fCrossSectionCalculator->Initialize(epsilon, E0_m, 0., fBeamPolarization,
                                      fTargetPolarization, 2);

  if(gamEnergy1 > lowestSecondaryEnergy)
  {
    // in interaction frame
    // calculate polarization transfer to the photon (in interaction plane)
    fFinalGammaPolarization = fCrossSectionCalculator->GetPol2();
    if(fVerboseLevel >= 1)
    {
      G4cout << " gammaPolarization1 = " << fFinalGammaPolarization << G4endl;
    }
    fFinalGammaPolarization.SetPhoton();

    // translate polarization into particle reference frame
    fFinalGammaPolarization.RotateAz(nInteractionFrame, gamDirection1);
    if(fFinalGammaPolarization.mag() > 1. + 1.e-8)
    {
      G4ExceptionDescription ed;
      ed << "ERROR in Polarizaed Compton Scattering !\n";
      ed << "Polarization of final photon more than 100%.\n";
      ed << fFinalGammaPolarization
         << " mag = " << fFinalGammaPolarization.mag() << "\n";
      G4Exception("G4PolarizedComptonModel::SampleSecondaries", "pol033",
                  FatalException, ed);
    }
    // store polarization vector
    fParticleChange->ProposePolarization(fFinalGammaPolarization);
    if(fVerboseLevel >= 1)
    {
      G4cout << " gammaPolarization1 = " << fFinalGammaPolarization << G4endl;
      G4cout << " GammaDirection1 = " << gamDirection1 << G4endl;
    }
  }

  // kinematic of the scattered electron
  G4double eKinEnergy = gamEnergy0 - gamEnergy1;

  if(eKinEnergy > lowestSecondaryEnergy)
  {
    G4ThreeVector eDirection =
      gamEnergy0 * gamDirection0 - gamEnergy1 * gamDirection1;
    eDirection = eDirection.unit();

    finalElectronPolarization = fCrossSectionCalculator->GetPol3();
    if(fVerboseLevel >= 1)
    {
      G4cout << " electronPolarization1 = " << finalElectronPolarization
             << G4endl;
    }
    // transfer into particle reference frame
    finalElectronPolarization.RotateAz(nInteractionFrame, eDirection);
    if(fVerboseLevel >= 1)
    {
      G4cout << " electronPolarization1 = " << finalElectronPolarization
             << G4endl << " ElecDirection = " << eDirection << G4endl;
    }

    // create G4DynamicParticle object for the electron.
    G4DynamicParticle* aElectron =
      new G4DynamicParticle(theElectron, eDirection, eKinEnergy);
    // store polarization vector
    if(finalElectronPolarization.mag() > 1. + 1.e-8)
    {
      G4ExceptionDescription ed;
      ed << "ERROR in Polarized Compton Scattering !\n";
      ed << "Polarization of final electron more than 100%.\n";
      ed << finalElectronPolarization
         << " mag = " << finalElectronPolarization.mag() << G4endl;
      G4Exception("G4PolarizedComptonModel::SampleSecondaries", "pol034",
                  FatalException, ed);
    }
    aElectron->SetPolarization(finalElectronPolarization.p1(),
                               finalElectronPolarization.p2(),
                               finalElectronPolarization.p3());
    fvect->push_back(aElectron);
  }
  else
  {
    edep += eKinEnergy;
  }
  // energy balance
  if(edep > 0.0)
  {
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedComptonModel::PrintWarning(const G4DynamicParticle* dp,
                                           G4int nloop, G4double grej,
                                           G4double onecos, G4double phi,
                                           const G4String sss) const
{
  G4ExceptionDescription ed;
  ed << "Problem of scattering sampling: " << sss << "\n"
     << "Niter= " << nloop << " grej= " << grej
     << " cos(theta)= " << 1.0 - onecos << " phi= " << phi << "\n"
     << "Gamma E(MeV)= " << dp->GetKineticEnergy() / MeV
     << " dir= " << dp->GetMomentumDirection()
     << " pol= " << dp->GetPolarization();
  G4Exception("G4PolarizedComptonModel::SampleSecondaries", "em0044",
              JustWarning, ed, "");
}

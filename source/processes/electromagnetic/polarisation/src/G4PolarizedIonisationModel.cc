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
// File name:     G4PolarizedIonisationModel
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Class Description:
//   Implementation of energy loss and delta-electron production by e+/e-
//   (including polarization effects)

#include "G4PolarizedIonisationModel.hh"

#include "G4ParticleChangeForLoss.hh"
#include "G4PhysicalConstants.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizedIonisationBhabhaXS.hh"
#include "G4PolarizedIonisationMollerXS.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedIonisationModel::G4PolarizedIonisationModel(
  const G4ParticleDefinition* p, const G4String& nam)
  : G4MollerBhabhaModel(p, nam)
  , fCrossSectionCalculator(nullptr)
{
  fBeamPolarization   = G4StokesVector::ZERO;
  fTargetPolarization = G4StokesVector::ZERO;

  fPositronPolarization = G4StokesVector::ZERO;
  fElectronPolarization = G4StokesVector::ZERO;

  isElectron = (p == theElectron);  // necessary due to wrong order in
                                    // G4MollerBhabhaModel constructor!

  if(!isElectron)
  {
    G4cout << " buildBhabha cross section " << isElectron << G4endl;
    fCrossSectionCalculator = new G4PolarizedIonisationBhabhaXS();
  }
  else
  {
    G4cout << " buildMoller cross section " << isElectron << G4endl;
    fCrossSectionCalculator = new G4PolarizedIonisationMollerXS();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedIonisationModel::~G4PolarizedIonisationModel()
{
  if(fCrossSectionCalculator)
  {
    delete fCrossSectionCalculator;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedIonisationModel::ComputeCrossSectionPerElectron(
  const G4ParticleDefinition* pd, G4double kinEnergy, G4double cut,
  G4double emax)
{
  G4double xs = G4MollerBhabhaModel::ComputeCrossSectionPerElectron(
    pd, kinEnergy, cut, emax);
  G4double factor = 1.;
  if(xs != 0.)
  {
    G4double tmax = MaxSecondaryEnergy(pd, kinEnergy);
    tmax          = std::min(emax, tmax);

    if(std::fabs(cut / emax - 1.) < 1.e-10)
      return xs;

    if(cut < tmax)
    {
      G4double xmin     = cut / kinEnergy;
      G4double xmax     = tmax / kinEnergy;
      G4double gam      = kinEnergy / electron_mass_c2 + 1.0;
      G4double crossPol = fCrossSectionCalculator->TotalXSection(
        xmin, xmax, gam, fBeamPolarization, fTargetPolarization);
      G4double crossUnpol = fCrossSectionCalculator->TotalXSection(
        xmin, xmax, gam, G4StokesVector::ZERO, G4StokesVector::ZERO);
      if(crossUnpol > 0.)
        factor = crossPol / crossUnpol;
    }
  }
  return xs * factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedIonisationModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* vdp, const G4MaterialCutsCouple*,
  const G4DynamicParticle* dp, G4double tmin, G4double maxEnergy)
{
  // *** obtain and save target and beam polarization ***
  G4PolarizationManager* polarizationManager =
    G4PolarizationManager::GetInstance();

  const G4Track* aTrack = fParticleChange->GetCurrentTrack();

  // obtain polarization of the beam
  fBeamPolarization = G4StokesVector(dp->GetPolarization());

  // obtain polarization of the media
  G4VPhysicalVolume* aPVolume    = aTrack->GetVolume();
  G4LogicalVolume* aLVolume      = aPVolume->GetLogicalVolume();
  const G4bool targetIsPolarized = polarizationManager->IsPolarized(aLVolume);
  fTargetPolarization = polarizationManager->GetVolumePolarization(aLVolume);

  // transfer target polarization in interaction frame
  if(targetIsPolarized)
    fTargetPolarization.rotateUz(dp->GetMomentumDirection());

  G4double tmax = std::min(maxEnergy, MaxSecondaryKinEnergy(dp));
  if(tmin >= tmax)
    return;

  G4double polL = fBeamPolarization.z() * fTargetPolarization.z();
  polL          = std::fabs(polL);
  G4double polT = fBeamPolarization.x() * fTargetPolarization.x() +
                  fBeamPolarization.y() * fTargetPolarization.y();
  polT = std::fabs(polT);

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double energy        = kineticEnergy + electron_mass_c2;
  G4double totalMomentum =
    std::sqrt(kineticEnergy * (energy + electron_mass_c2));
  G4double xmin   = tmin / kineticEnergy;
  G4double xmax   = tmax / kineticEnergy;
  G4double gam    = energy / electron_mass_c2;
  G4double gamma2 = gam * gam;
  G4double gmo    = gam - 1.;
  G4double gmo2   = gmo * gmo;
  G4double gmo3   = gmo2 * gmo;
  G4double gpo    = gam + 1.;
  G4double gpo2   = gpo * gpo;
  G4double gpo3   = gpo2 * gpo;
  G4double x, y, q, grej, grej2;
  G4double z  = 0.;
  G4double xs = 0., phi = 0.;
  G4ThreeVector direction = dp->GetMomentumDirection();

  //(Polarized) Moller (e-e-) scattering
  if(isElectron)
  {
    // *** dice according to polarized cross section
    G4double G = ((2.0 * gam - 1.0) / gamma2) * (1. - polT - polL * gam);
    G4double H = (sqr(gam - 1.0) / gamma2) *
                 (1. + polT + polL * ((gam + 3.) / (gam - 1.)));

    y     = 1.0 - xmax;
    grej  = 1.0 - G * xmax + xmax * xmax * (H + (1.0 - G * y) / (y * y));
    grej2 = 1.0 - G * xmin + xmin * xmin * (H + (1.0 - G * y) / (y * y));
    if(grej2 > grej)
      grej = grej2;
    G4double prefM = gamma2 * classic_electr_radius * classic_electr_radius /
                     (gmo2 * (gam + 1.0));
    grej *= prefM;
    do
    {
      q = G4UniformRand();
      x = xmin * xmax / (xmin * (1.0 - q) + xmax * q);
      if(fCrossSectionCalculator)
      {
        fCrossSectionCalculator->Initialize(x, gam, phi, fBeamPolarization,
                                            fTargetPolarization, 1);
        xs = fCrossSectionCalculator->XSection(G4StokesVector::ZERO,
                                               G4StokesVector::ZERO);
        z  = xs * sqr(x) * 4.;
        if(grej < z)
        {
          G4ExceptionDescription ed;
          ed << "WARNING : error in Moller rejection routine! \n"
             << " z = " << z << " grej=" << grej << "\n";
          G4Exception("G4PolarizedIonisationModel::SampleSecondaries", "pol019",
                      JustWarning, ed);
        }
      }
      else
      {
        G4cout << "No calculator in Moller scattering" << G4endl;
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(grej * G4UniformRand() > z);
    // Bhabha (e+e-) scattering
  }
  else
  {
    // *** dice according to polarized cross section
    y    = xmax * xmax;
    grej = 0.;
    grej += y * y * gmo3 * (1. + (polL + polT) * (gam + 3.) / gmo);
    grej += -2. * xmin * xmin * xmin * gam * gmo2 *
            (1. - (polL + polT) * (gam + 3.) / gmo);
    grej += y * y * gmo * (3. * gamma2 + 6. * gam + 4.) *
            (1. + (polL * (3. * gam + 1.) * (gamma2 + gam + 1.) +
                   polT * ((gam + 2.) * gamma2 + 1.)) /
                    (gmo * (3. * gam * (gam + 2.) + 4.)));
    grej /= gpo3;
    grej += -xmin * (2. * gamma2 + 4. * gam + 1.) *
            (1. - gam * (polL * (2. * gam + 1.) + polT) /
                    (2. * gam * (gam + 2.) + 1.)) /
            gpo2;
    grej += gamma2 / (gamma2 - 1.);
    G4double prefB =
      classic_electr_radius * classic_electr_radius / (gam - 1.0);
    grej *= prefB;

    do
    {
      q = G4UniformRand();
      x = xmin * xmax / (xmin * (1.0 - q) + xmax * q);
      if(fCrossSectionCalculator)
      {
        fCrossSectionCalculator->Initialize(x, gam, phi, fBeamPolarization,
                                            fTargetPolarization, 1);
        xs = fCrossSectionCalculator->XSection(G4StokesVector::ZERO,
                                               G4StokesVector::ZERO);
        z  = xs * sqr(x) * 4.;
      }
      else
      {
        G4cout << "No calculator in Bhabha scattering" << G4endl;
      }

      if(z > grej)
      {
        G4ExceptionDescription ed;
        ed << "G4PolarizedIonisationModel::SampleSecondaries Warning!\n "
           << "Majorant " << grej << " < " << z << " for x= " << x << G4endl
           << " e+e- (Bhabha) scattering"
           << " at KinEnergy " << kineticEnergy << G4endl;
        G4Exception("G4PolarizedIonisationModel::SampleSecondaries", "pol020",
                    JustWarning, ed);
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(grej * G4UniformRand() > z);
  }

  // polar asymmetries (due to transverse polarizations)
  if(fCrossSectionCalculator)
  {
    grej = xs * 2.;
    do
    {
      phi = twopi * G4UniformRand();
      fCrossSectionCalculator->Initialize(x, gam, phi, fBeamPolarization,
                                          fTargetPolarization, 1);
      xs = fCrossSectionCalculator->XSection(G4StokesVector::ZERO,
                                             G4StokesVector::ZERO);
      if(xs > grej)
      {
        if(isElectron)
        {
          G4ExceptionDescription ed;
          ed << "Majorant " << grej << " < " << xs << " for phi= " << phi
             << "\n"
             << " e-e- (Moller) scattering\n"
             << "PHI DICING\n";
          G4Exception("G4PolarizedIonisationModel::SampleSecondaries", "pol021",
                      JustWarning, ed);
        }
        else
        {
          G4ExceptionDescription ed;
          ed << "Majorant " << grej << " < " << xs << " for phi= " << phi
             << "\n"
             << " e+e- (Bhabha) scattering\n"
             << "PHI DICING\n";
          G4Exception("G4PolarizedIonisationModel::SampleSecondaries", "pol022",
                      JustWarning, ed);
        }
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(grej * G4UniformRand() > xs);
  }

  // fix kinematics of delta electron
  G4double deltaKinEnergy = x * kineticEnergy;
  G4double deltaMomentum =
    std::sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0 * electron_mass_c2));
  G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
                  (deltaMomentum * totalMomentum);
  G4double sint = 1.0 - cost * cost;
  if(sint > 0.0)
    sint = std::sqrt(sint);

  G4ThreeVector deltaDirection(-sint * std::cos(phi), -sint * std::sin(phi),
                               cost);
  deltaDirection.rotateUz(direction);

  // primary change
  kineticEnergy -= deltaKinEnergy;
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);

  if(kineticEnergy > DBL_MIN)
  {
    G4ThreeVector dir =
      totalMomentum * direction - deltaMomentum * deltaDirection;
    direction = dir.unit();
    fParticleChange->SetProposedMomentumDirection(direction);
  }

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta =
    new G4DynamicParticle(theElectron, deltaDirection, deltaKinEnergy);
  vdp->push_back(delta);

  // get interaction frame
  G4ThreeVector nInteractionFrame =
    G4PolarizationHelper::GetFrame(direction, deltaDirection);

  if(fCrossSectionCalculator)
  {
    // calculate mean final state polarizations
    fBeamPolarization.InvRotateAz(nInteractionFrame, direction);
    fTargetPolarization.InvRotateAz(nInteractionFrame, direction);
    fCrossSectionCalculator->Initialize(x, gam, phi, fBeamPolarization,
                                        fTargetPolarization, 2);

    // electron/positron
    fPositronPolarization = fCrossSectionCalculator->GetPol2();
    fPositronPolarization.RotateAz(nInteractionFrame, direction);

    fParticleChange->ProposePolarization(fPositronPolarization);

    // electron
    fElectronPolarization = fCrossSectionCalculator->GetPol3();
    fElectronPolarization.RotateAz(nInteractionFrame, deltaDirection);
    delta->SetPolarization(fElectronPolarization.x(), fElectronPolarization.y(),
                           fElectronPolarization.z());
  }
  else
  {
    fPositronPolarization = G4StokesVector::ZERO;
    fElectronPolarization = G4StokesVector::ZERO;
  }
}

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
// File name:     G4PolarizedAnnihilationModel
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   Implementation of polarized gamma Annihilation scattering on free electron

#include "G4PolarizedAnnihilationModel.hh"

#include "G4Gamma.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PhysicalConstants.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizedAnnihilationXS.hh"
#include "G4StokesVector.hh"
#include "G4TrackStatus.hh"

G4PolarizedAnnihilationModel::G4PolarizedAnnihilationModel(
  const G4ParticleDefinition* p, const G4String& nam)
  : G4eeToTwoGammaModel(p, nam)
  , fCrossSectionCalculator(nullptr)
  , fParticleChange(nullptr)
  , fVerboseLevel(0)
{
  fCrossSectionCalculator  = new G4PolarizedAnnihilationXS();
  fBeamPolarization        = G4StokesVector::ZERO;
  fTargetPolarization      = G4StokesVector::ZERO;
  fFinalGamma1Polarization = G4StokesVector::ZERO;
  fFinalGamma2Polarization = G4StokesVector::ZERO;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedAnnihilationModel::~G4PolarizedAnnihilationModel()
{
  delete fCrossSectionCalculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedAnnihilationModel::Initialise(const G4ParticleDefinition* part,
                                              const G4DataVector& dv)
{
  G4eeToTwoGammaModel::Initialise(part, dv);
  if(fParticleChange)
  {
    return;
  }
  fParticleChange = GetParticleChangeForGamma();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedAnnihilationModel::ComputeCrossSectionPerElectron(
  G4double kinEnergy)
{
  // cross section from base model
  G4double xs = G4eeToTwoGammaModel::ComputeCrossSectionPerElectron(kinEnergy);

  G4double polzz = fBeamPolarization.z() * fTargetPolarization.z();
  G4double poltt = fBeamPolarization.x() * fTargetPolarization.x() +
                   fBeamPolarization.y() * fTargetPolarization.y();
  if(polzz != 0 || poltt != 0)
  {
    G4double xval, lasym, tasym;
    ComputeAsymmetriesPerElectron(kinEnergy, xval, lasym, tasym);
    xs *= (1. + polzz * lasym + poltt * tasym);
  }

  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedAnnihilationModel::ComputeAsymmetriesPerElectron(
  G4double ene, G4double& valueX, G4double& valueA, G4double& valueT)
{
  // *** calculate asymmetries
  G4double gam = 1. + ene / electron_mass_c2;
  G4double xs0 = fCrossSectionCalculator->TotalXSection(
    0., 1., gam, G4StokesVector::ZERO, G4StokesVector::ZERO);
  G4double xsA = fCrossSectionCalculator->TotalXSection(
    0., 1., gam, G4StokesVector::P3, G4StokesVector::P3);
  G4double xsT1 = fCrossSectionCalculator->TotalXSection(
    0., 1., gam, G4StokesVector::P1, G4StokesVector::P1);
  G4double xsT2 = fCrossSectionCalculator->TotalXSection(
    0., 1., gam, G4StokesVector::P2, G4StokesVector::P2);
  G4double xsT = 0.5 * (xsT1 + xsT2);

  valueX = xs0;
  valueA = xsA / xs0 - 1.;
  valueT = xsT / xs0 - 1.;

  if((valueA < -1) || (1 < valueA))
  {
    G4ExceptionDescription ed;
    ed << " ERROR PolarizedAnnihilationPS::ComputeAsymmetries \n";
    ed << " something wrong in total cross section calculation (valueA)\n";
    ed << " LONG: " << valueX << "\t" << valueA << "\t" << valueT
       << "   energy = " << gam << G4endl;
    G4Exception("G4PolarizedAnnihilationModel::ComputeAsymmetriesPerElectron",
                "pol004", JustWarning, ed);
  }
  if((valueT < -1) || (1 < valueT))
  {
    G4ExceptionDescription ed;
    ed << " ERROR PolarizedAnnihilationPS::ComputeAsymmetries \n";
    ed << " something wrong in total cross section calculation (valueT)\n";
    ed << " TRAN: " << valueX << "\t" << valueA << "\t" << valueT
       << "   energy = " << gam << G4endl;
    G4Exception("G4PolarizedAnnihilationModel::ComputeAsymmetriesPerElectron",
                "pol005", JustWarning, ed);
  }
}

void G4PolarizedAnnihilationModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* fvect, const G4MaterialCutsCouple*,
  const G4DynamicParticle* dp, G4double, G4double)
{
  const G4Track* aTrack = fParticleChange->GetCurrentTrack();

  // kill primary
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);

  // V.Ivanchenko add protection against zero kin energy
  G4double PositKinEnergy = dp->GetKineticEnergy();

  if(PositKinEnergy == 0.0)
  {
    G4double cosTeta = 2. * G4UniformRand() - 1.;
    G4double sinTeta = std::sqrt((1.0 - cosTeta) * (1.0 + cosTeta));
    G4double phi     = twopi * G4UniformRand();
    G4ThreeVector dir(sinTeta * std::cos(phi), sinTeta * std::sin(phi),
                      cosTeta);
    fvect->push_back(
      new G4DynamicParticle(G4Gamma::Gamma(), dir, electron_mass_c2));
    fvect->push_back(
      new G4DynamicParticle(G4Gamma::Gamma(), -dir, electron_mass_c2));
    return;
  }

  // *** obtain and save target and beam polarization ***
  G4PolarizationManager* polarizationManager =
    G4PolarizationManager::GetInstance();

  // obtain polarization of the beam
  fBeamPolarization = G4StokesVector(aTrack->GetPolarization());

  // obtain polarization of the media
  G4VPhysicalVolume* aPVolume    = aTrack->GetVolume();
  G4LogicalVolume* aLVolume      = aPVolume->GetLogicalVolume();
  const G4bool targetIsPolarized = polarizationManager->IsPolarized(aLVolume);
  fTargetPolarization = polarizationManager->GetVolumePolarization(aLVolume);

  if(fVerboseLevel >= 1)
  {
    G4cout << "G4PolarizedComptonModel::SampleSecondaries in "
           << aLVolume->GetName() << G4endl;
  }

  // transfer target electron polarization in frame of positron
  if(targetIsPolarized)
    fTargetPolarization.rotateUz(dp->GetMomentumDirection());

  G4ParticleMomentum PositDirection = dp->GetMomentumDirection();

  // polar asymmetry:
  G4double polarization = fBeamPolarization.p3() * fTargetPolarization.p3();

  G4double gamam1 = PositKinEnergy / electron_mass_c2;
  G4double gama = gamam1 + 1., gamap1 = gamam1 + 2.;
  G4double sqgrate = std::sqrt(gamam1 / gamap1) / 2.,
           sqg2m1  = std::sqrt(gamam1 * gamap1);

  // limits of the energy sampling
  G4double epsilmin = 0.5 - sqgrate, epsilmax = 0.5 + sqgrate;
  G4double epsilqot = epsilmax / epsilmin;

  // sample the energy rate of the created gammas
  // note: for polarized partices, the actual dicing strategy
  //       will depend on the energy, and the degree of polarization !!
  G4double epsil;
  G4double gmax = 1. + std::fabs(polarization);  // crude estimate

  fCrossSectionCalculator->Initialize(epsilmin, gama, 0., fBeamPolarization,
                                      fTargetPolarization);
  if(fCrossSectionCalculator->DiceEpsilon() < 0.)
  {
    G4ExceptionDescription ed;
    ed << "ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
       << "epsilmin DiceRoutine not appropriate ! "
       << fCrossSectionCalculator->DiceEpsilon() << G4endl;
    G4Exception("G4PolarizedAnnihilationModel::SampleSecondaries", "pol006",
                JustWarning, ed);
  }

  fCrossSectionCalculator->Initialize(epsilmax, gama, 0., fBeamPolarization,
                                      fTargetPolarization);
  if(fCrossSectionCalculator->DiceEpsilon() < 0)
  {
    G4ExceptionDescription ed;
    ed << "ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
       << "epsilmax DiceRoutine not appropriate ! "
       << fCrossSectionCalculator->DiceEpsilon() << G4endl;
    G4Exception("G4PolarizedAnnihilationModel::SampleSecondaries", "pol007",
                JustWarning, ed);
  }

  G4int ncount        = 0;
  G4double trejectmax = 0.;
  G4double treject;

  do
  {
    epsil = epsilmin * std::pow(epsilqot, G4UniformRand());

    fCrossSectionCalculator->Initialize(epsil, gama, 0., fBeamPolarization,
                                        fTargetPolarization, 1);

    treject = fCrossSectionCalculator->DiceEpsilon();
    treject *= epsil;

    if(treject > gmax || treject < 0.)
    {
      G4ExceptionDescription ed;
      ed << "ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
         << " eps (" << epsil
         << ") rejection does not work properly: " << treject << G4endl;
      G4Exception("G4PolarizedAnnihilationModel::SampleSecondaries", "pol008",
                  JustWarning, ed);
    }
    ++ncount;
    if(treject > trejectmax)
      trejectmax = treject;
    if(ncount > 1000)
    {
      G4ExceptionDescription ed;
      ed << "WARNING  in PolarizedAnnihilationPS::PostStepDoIt\n"
         << "eps dicing very inefficient =" << trejectmax / gmax << ", "
         << treject / gmax << ".  For secondary energy = " << epsil << "   "
         << ncount << G4endl;
      G4Exception("G4PolarizedAnnihilationModel::SampleSecondaries", "pol009",
                  JustWarning, ed);
      break;
    }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while(treject < gmax * G4UniformRand());

  // scattered Gamma angles. ( Z - axis along the parent positron)
  G4double cost = (epsil * gamap1 - 1.) / (epsil * sqg2m1);
  G4double sint = std::sqrt((1. + cost) * (1. - cost));
  G4double phi  = 0.;
  G4double beamTrans =
    std::sqrt(sqr(fBeamPolarization.p1()) + sqr(fBeamPolarization.p2()));
  G4double targetTrans =
    std::sqrt(sqr(fTargetPolarization.p1()) + sqr(fTargetPolarization.p2()));

  do
  {
    phi = twopi * G4UniformRand();
    fCrossSectionCalculator->Initialize(epsil, gama, 0., fBeamPolarization,
                                        fTargetPolarization, 2);

    G4double gdiced = fCrossSectionCalculator->getVar(0);
    gdiced += fCrossSectionCalculator->getVar(3) * fBeamPolarization.p3() *
              fTargetPolarization.p3();
    gdiced += 1. *
              (std::fabs(fCrossSectionCalculator->getVar(1)) +
               std::fabs(fCrossSectionCalculator->getVar(2))) *
              beamTrans * targetTrans;
    gdiced += 1. * std::fabs(fCrossSectionCalculator->getVar(4)) *
              (std::fabs(fBeamPolarization.p3()) * targetTrans +
               std::fabs(fTargetPolarization.p3()) * beamTrans);

    G4double gdist = fCrossSectionCalculator->getVar(0);
    gdist += fCrossSectionCalculator->getVar(3) * fBeamPolarization.p3() *
             fTargetPolarization.p3();
    gdist += fCrossSectionCalculator->getVar(1) *
             (std::cos(phi) * fBeamPolarization.p1() +
              std::sin(phi) * fBeamPolarization.p2()) *
             (std::cos(phi) * fTargetPolarization.p1() +
              std::sin(phi) * fTargetPolarization.p2());
    gdist += fCrossSectionCalculator->getVar(2) *
             (std::cos(phi) * fBeamPolarization.p2() -
              std::sin(phi) * fBeamPolarization.p1()) *
             (std::cos(phi) * fTargetPolarization.p2() -
              std::sin(phi) * fTargetPolarization.p1());
    gdist +=
      fCrossSectionCalculator->getVar(4) *
      (std::cos(phi) * fBeamPolarization.p3() * fTargetPolarization.p1() +
       std::cos(phi) * fBeamPolarization.p1() * fTargetPolarization.p3() +
       std::sin(phi) * fBeamPolarization.p3() * fTargetPolarization.p2() +
       std::sin(phi) * fBeamPolarization.p2() * fTargetPolarization.p3());

    treject = gdist / gdiced;
    if(treject > 1. + 1.e-10 || treject < 0)
    {
      G4ExceptionDescription ed;
      ed << "!!!ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
         << " phi rejection does not work properly: " << treject << G4endl;
      G4cout << " gdiced = " << gdiced << G4endl;
      G4cout << " gdist = " << gdist << G4endl;
      G4cout << " epsil = " << epsil << G4endl;
      G4Exception("G4PolarizedAnnihilationModel::SampleSecondaries", "pol009",
                  JustWarning, ed);
    }

    if(treject < 1.e-3)
    {
      G4ExceptionDescription ed;
      ed << "!!!ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
         << " phi rejection does not work properly: " << treject << "\n";
      G4cout << " gdiced=" << gdiced << "   gdist=" << gdist << "\n";
      G4cout << " epsil = " << epsil << G4endl;
      G4Exception("G4PolarizedAnnihilationModel::SampleSecondaries", "pol010",
                  JustWarning, ed);
    }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while(treject < G4UniformRand());

  G4double dirx = sint * std::cos(phi);
  G4double diry = sint * std::sin(phi);
  G4double dirz = cost;

  // kinematic of the created pair
  G4double TotalAvailableEnergy = PositKinEnergy + 2 * electron_mass_c2;
  G4double Phot1Energy          = epsil * TotalAvailableEnergy;
  G4double Phot2Energy          = (1. - epsil) * TotalAvailableEnergy;

  // *** prepare calculation of polarization transfer ***
  G4ThreeVector Phot1Direction(dirx, diry, dirz);

  // get interaction frame
  G4ThreeVector nInteractionFrame =
    G4PolarizationHelper::GetFrame(PositDirection, Phot1Direction);

  // define proper in-plane and out-of-plane component of initial spins
  fBeamPolarization.InvRotateAz(nInteractionFrame, PositDirection);
  fTargetPolarization.InvRotateAz(nInteractionFrame, PositDirection);

  // calculate spin transfere matrix

  fCrossSectionCalculator->Initialize(epsil, gama, phi, fBeamPolarization,
                                      fTargetPolarization, 2);

  Phot1Direction.rotateUz(PositDirection);
  // create G4DynamicParticle object for the particle1
  G4DynamicParticle* aParticle1 =
    new G4DynamicParticle(G4Gamma::Gamma(), Phot1Direction, Phot1Energy);
  fFinalGamma1Polarization = fCrossSectionCalculator->GetPol2();
  G4double n1              = fFinalGamma1Polarization.mag2();
  if(n1 > 1.)
  {
    G4ExceptionDescription ed;
    ed << "ERROR: PolarizedAnnihilation Polarization Vector at epsil = "
       << epsil << " is too large!!! \n"
       << "annihi pol1= " << fFinalGamma1Polarization << ", (" << n1 << ")\n";
    fFinalGamma1Polarization *= 1. / std::sqrt(n1);
    G4Exception("G4PolarizedAnnihilationModel::SampleSecondaries", "pol011",
                JustWarning, ed);
  }

  // define polarization of first final state photon
  fFinalGamma1Polarization.SetPhoton();
  fFinalGamma1Polarization.RotateAz(nInteractionFrame, Phot1Direction);
  aParticle1->SetPolarization(fFinalGamma1Polarization.p1(),
                              fFinalGamma1Polarization.p2(),
                              fFinalGamma1Polarization.p3());

  fvect->push_back(aParticle1);

  // **********************************************************************

  G4double Eratio = Phot1Energy / Phot2Energy;
  G4double PositP =
    std::sqrt(PositKinEnergy * (PositKinEnergy + 2. * electron_mass_c2));
  G4ThreeVector Phot2Direction(-dirx * Eratio, -diry * Eratio,
                               (PositP - dirz * Phot1Energy) / Phot2Energy);
  Phot2Direction.rotateUz(PositDirection);
  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2 =
    new G4DynamicParticle(G4Gamma::Gamma(), Phot2Direction, Phot2Energy);

  // define polarization of second final state photon
  fFinalGamma2Polarization = fCrossSectionCalculator->GetPol3();
  G4double n2              = fFinalGamma2Polarization.mag2();
  if(n2 > 1.)
  {
    G4ExceptionDescription ed;
    ed << "ERROR: PolarizedAnnihilation Polarization Vector at epsil = "
       << epsil << " is too large!!! \n";
    ed << "annihi pol2= " << fFinalGamma2Polarization << ", (" << n2 << ")\n";

    G4Exception("G4PolarizedAnnihilationModel::SampleSecondaries", "pol012",
                JustWarning, ed);
    fFinalGamma2Polarization *= 1. / std::sqrt(n2);
  }
  fFinalGamma2Polarization.SetPhoton();
  fFinalGamma2Polarization.RotateAz(nInteractionFrame, Phot2Direction);
  aParticle2->SetPolarization(fFinalGamma2Polarization.p1(),
                              fFinalGamma2Polarization.p2(),
                              fFinalGamma2Polarization.p3());

  fvect->push_back(aParticle2);
}

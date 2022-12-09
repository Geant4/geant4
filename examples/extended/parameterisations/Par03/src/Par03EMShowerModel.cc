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
#include "Par03EMShowerModel.hh"
#include "Par03EMShowerMessenger.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4FastHit.hh"
#include "Randomize.hh"
#include "G4FastSimHitMaker.hh"

Par03EMShowerModel::Par03EMShowerModel(G4String aModelName, G4Region* aEnvelope)
  : G4VFastSimulationModel(aModelName, aEnvelope)
  , fMessenger(new Par03EMShowerMessenger(this))
  , fHitMaker(new G4FastSimHitMaker)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03EMShowerModel::Par03EMShowerModel(G4String aModelName)
  : G4VFastSimulationModel(aModelName)
  , fMessenger(new Par03EMShowerMessenger(this))
  , fHitMaker(new G4FastSimHitMaker)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03EMShowerModel::~Par03EMShowerModel() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par03EMShowerModel::IsApplicable(
  const G4ParticleDefinition& aParticleType)
{
  return &aParticleType == G4Electron::ElectronDefinition() ||
         &aParticleType == G4Positron::PositronDefinition() ||
         &aParticleType == G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par03EMShowerModel::ModelTrigger(const G4FastTrack& aFastTrack)
{
  // Check energy
  if(aFastTrack.GetPrimaryTrack()->GetKineticEnergy() < 1 * GeV)
  {
    return false;
  }
  // Check length of detector
  // Calculate depth of the detector along shower axis to verify if shower
  // will fit inside. Required max shower depth is defined by fLongMaxDepth, and
  // can be changed with UI command `/Par03/fastSim/longitudinalProfile/maxDepth
  G4double X0 = aFastTrack.GetPrimaryTrack()->GetMaterial()->GetRadlen();
  auto particleDirection     = aFastTrack.GetPrimaryTrackLocalDirection();
  auto particlePosition      = aFastTrack.GetPrimaryTrackLocalPosition();
  G4double detectorDepthInMM = aFastTrack.GetEnvelopeSolid()->DistanceToOut(
    particlePosition, particleDirection);
  G4double detectorDepthInX0 = detectorDepthInMM / X0;
  // check if detector depth is sufficient to create showers
  if(detectorDepthInX0 < fLongMaxDepth)
  {
    return false;
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03EMShowerModel::DoIt(const G4FastTrack& aFastTrack,
                              G4FastStep& aFastStep)
{
  // Remove particle from further processing by G4
  aFastStep.KillPrimaryTrack();
  aFastStep.SetPrimaryTrackPathLength(0.0);
  G4double energy = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();
  // No need to create any deposit, it will be handled by this model (and
  // G4FastSimHitMaker that will call the sensitive detector)
  aFastStep.SetTotalEnergyDeposited(0);
  auto particlePosition  = aFastTrack.GetPrimaryTrackLocalPosition();
  auto particleDirection = aFastTrack.GetPrimaryTrackLocalDirection();

  // Calculate how to create energy deposits
  // Following PDG 33.5 chapter
  // material calculation assumes homogeneous detector (true for Par03 example)
  auto material       = aFastTrack.GetPrimaryTrack()->GetMaterial();
  G4double materialX0 = material->GetRadlen();
  G4double materialZ  = material->GetZ();
  // EC estimation follows PDG fit to solids in Fig. 33.14 (rms 2.2%)
  G4double materialEc = 610 * MeV / (materialZ + 1.24);
  // RM estimation follows PDG Eq. (33.37) (rms 2.2%)
  G4double materialRM = 21.2052 * MeV * materialX0 / materialEc;
  G4double particleY  = energy / materialEc;
  // Estimate shower maximum and alpha parameter of Gamma distribution
  // that describes the longitudinal profile (PDG Eq. (33.35))
  // unless alpha is specified by UI command
  if(fAlpha < 0)
  {
    // from PDG Eq. (33.36)
    G4double particleTmax = std::log(particleY);
    if(aFastTrack.GetPrimaryTrack()->GetParticleDefinition() ==
       G4Gamma::GammaDefinition())
    {
      particleTmax += 0.5;
    }
    else
    {
      particleTmax -= 0.5;
    }
    fAlpha = particleTmax * fBeta + 1;
  }
  // Unless sigma of Gaussian distribution describing the transverse profile
  // is specified by UI command, use value calculated from Moliere Radius
  if(fSigma < 0)
  {
    // 90% of shower is contained within 1 * R_M
    // 1.645 * std dev of Gaussian contains 90%
    fSigma = materialRM / 1.645;
  }

  // Calculate rotation matrix along the particle momentum direction
  // It will rotate the shower axes to match the incoming particle direction
  G4RotationMatrix rotMatrix = G4RotationMatrix();
  double particleTheta       = particleDirection.theta();
  double particlePhi         = particleDirection.phi();
  double epsilon             = 1e-3;
  rotMatrix.rotateY(particleTheta);
  // do not use (random) phi if x==y==0
  if(!(std::fabs(particleDirection.x()) < epsilon &&
       std::fabs(particleDirection.y()) < epsilon))
    rotMatrix.rotateZ(particlePhi);

  // Create hits
  // First use rejecton sampling to sample from Gamma distribution
  // then get random numbers from uniform distribution for azimuthal angle, and
  // from Gaussian for radius
  G4ThreeVector position;
  G4double gammaMax   = Gamma((fAlpha - 1) / fBeta, fAlpha, fBeta);
  G4int generatedHits = 0;
  while(generatedHits < fNbOfHits)
  {
    G4double random1 = G4UniformRand() * fLongMaxDepth;
    G4double random2 = G4UniformRand() * gammaMax;
    if(Gamma(random1, fAlpha, fBeta) >= random2)
    {
      // Generate corresponding rho (phi) from Gaussian (flat) distribution
      G4double phiPosition = G4UniformRand() * 2 * CLHEP::pi;
      G4double rhoPosition = G4RandGauss::shoot(0, fSigma);
      position             = particlePosition +
                 rotMatrix * G4ThreeVector(rhoPosition * std::sin(phiPosition),
                                           rhoPosition * std::cos(phiPosition),
                                           random1 * materialX0);
      // Create energy deposit in the detector
      // This will call appropriate sensitive detector class
      fHitMaker->make(G4FastHit(position, energy / fNbOfHits), aFastTrack);
      generatedHits++;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03EMShowerModel::Print() const
{
  G4cout << "Par03EMShowerModel: " << G4endl;
  G4cout << "Gaussian distribution (transverse plane): \tmu = 0, sigma = "
         << G4BestUnit(fSigma, "Length") << G4endl;
  if(fSigma < 0)
    G4cout << "Negative sigma value means that it will be recalculated "
              "from the value of the Moliere radius of the detector material, "
              "taking into account that 90% of the area below the Gaussian "
              "distribution (from mu - 1.645 sigma to mu + 1.645 sigma) "
              "corresponds to area within 1 Moliere radius."
           << G4endl;
  G4cout << "Gamma distribution (along shower axis): \talpha = " << fAlpha
         << ", beta = " << fBeta << ", max depth = " << fLongMaxDepth << " X0"
         << G4endl;
  if(fAlpha < 0)
    G4cout << "Negative alpha value means that it will be recalculated "
              "from the critical energy of the detector material, particle "
              "type, and beta parameter.\n alpha = beta * T_max, where T_max = "
              "ln(E/E_C) + C\n where E is particle energy, E_C is critical "
              "energy in the material, and constant C = -0.5 for electrons and "
              "0.5 for photons (Eq. (33.36) from PDG)."
           << G4endl;
  G4cout << "Number of created energy deposits: " << fNbOfHits << G4endl;
}
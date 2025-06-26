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

#include "G4ChemReboundTransportation.hh"

#include "G4DNAMolecularMaterial.hh"
#include "G4H3O.hh"
#include "G4ITNavigator.hh"
#include "G4ITSafetyHelper.hh"  // Not used yet
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4Molecule.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4RandomDirection.hh"
#include "G4SafetyHelper.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackingInformation.hh"
#include "G4UnitsTable.hh"
#include "G4VUserBrownianAction.hh"
#include "Randomize.hh"

#include <CLHEP/Random/Stat.h>
using namespace std;

#ifndef State
#  define State(theXInfo) (GetState<G4ITBrownianState>()->theXInfo)
#endif

static G4double InvErfc(G4double x)
{
  return CLHEP::HepStat::inverseErf(1. - x);
}

#ifndef State
#  define State(theXInfo) (GetState<G4ITTransportationState>()->theXInfo)
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChemReboundTransportation::G4ChemReboundTransportation(const G4String& aName,
                                                         const G4DNABoundingBox* pB,
                                                         G4int verbosity)
  : G4ITTransportation(aName, verbosity), fpBoundingBox(pB)
{
  fVerboseLevel = 0;
  fpState = std::make_shared<G4ITBrownianState>();
  SetProcessSubType(fLowEnergyBrownianTransportation);
  fNistWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  fInternalMinTimeStep = 1 * CLHEP::ps;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChemReboundTransportation::G4ITBrownianState::G4ITBrownianState()
{
  fTimeStepReachedLimit = false;
  fRandomNumber = -1;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChemReboundTransportation::StartTracking(G4Track* track)
{
  fpState = std::make_shared<G4ITBrownianState>();
  SetInstantiateProcessState(false);
  G4ITTransportation::StartTracking(track);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChemReboundTransportation::BuildPhysicsTable(const G4ParticleDefinition& particle)
{
  fpSafetyHelper->InitialiseHelper();
  G4ITTransportation::BuildPhysicsTable(particle);

  if (fpBoundingBox == nullptr) {
    G4ExceptionDescription errMsg;
    errMsg << "fpBoundingBox is nullptr";
    G4Exception(
      "ChemReboundTransportation::BuildPhysicsTable"
      "ChemReboundTransportation",
      "ChemReboundTransportation", FatalErrorInArgument, errMsg);
  }
  G4double halfSize =
    std::min({fpBoundingBox->halfSideLengthInX(), fpBoundingBox->halfSideLengthInY(),
              fpBoundingBox->halfSideLengthInZ()});
  fMaximumTimeStep = (halfSize * halfSize) / (60 * G4H3O::Definition()->GetDiffusionCoefficient());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChemReboundTransportation::ComputeStep(const G4Track& track, const G4Step& step,
                                              const G4double timeStep, G4double& spaceStep)
{
  if (GetIT(track)->GetTrackingInfo()->IsLeadingStep()) {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "ComputeStep is called while the track has"
                            "the minimum interaction time";
    exceptionDescription << " so it should not recompute a timeStep ";
    G4Exception("ChemReboundTransportation::ComputeStep", "ChemReboundTransportation0001",
                FatalErrorInArgument, exceptionDescription);
  }
  State(fGeometryLimitedStep) = false;
  if (timeStep == 0) {
    State(fTransportEndPosition) = track.GetPosition();
    spaceStep = 0.;
  }
  else {
    auto molConf = GetMolecule(track)->GetMolecularConfiguration();
    spaceStep = calculateDistanceFromTimeStep(molConf, timeStep);
  }
  State(fTransportEndPosition) =
    BouncingAction(track.GetPosition() + spaceStep * G4RandomDirection());
  State(fEndPointDistance) = (track.GetPosition() - State(fTransportEndPosition)).mag();
  if (fVerboseLevel > 1)
  // if(GetMolecule(track)->GetName() == "e_aq^-1")
  {
    G4cout << G4endl;
    G4cout << "ComputeStep: timeStep : " << G4BestUnit(timeStep, "Time")
           << "  State(theInteractionTimeLeft) : " << State(theInteractionTimeLeft)
           << "  State(fEndPointDistance) : " << G4BestUnit(State(fEndPointDistance), "Length")
           << " trackID : " << track.GetTrackID()
           << " Molecule name: " << "track.GetPosition() : " << track.GetPosition()
           << "  State(fTransportEndPosition) : " << State(fTransportEndPosition) << "  "
           << GetMolecule(track)->GetName() << "   Diffusion length : " << G4endl;
  }
  State(fCandidateEndGlobalTime) = step.GetPreStepPoint()->GetGlobalTime() + timeStep;
  State(fEndGlobalTimeComputed) = true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4ChemReboundTransportation::PostStepDoIt(const G4Track& track,
                                                             const G4Step& step)
{
  G4ITTransportation::PostStepDoIt(track, step);

#ifdef G4VERBOSE
  //    DEBUG
  if (fVerboseLevel > 1)
    // if(GetMolecule(track)->GetName() == "e_aq^-1")
    if (GetIT(track)->GetTrackingInfo()->IsLeadingStep()) {
      G4cout << "ChemReboundTransportation::PostStepDoIt() :" << " trackID : " << track.GetTrackID()
             << " Molecule name: " << "prePosition : " << step.GetPreStepPoint()->GetPosition()
             << "  postPostion : " << step.GetPostStepPoint()->GetPosition() << "  "
             << GetMolecule(track)->GetName()
             << "   Diffusion length : " << G4BestUnit(step.GetStepLength(), "Length")
             << " within time step : " << G4BestUnit(step.GetDeltaTime(), "Time")
             << "\t Current global time : " << G4BestUnit(track.GetGlobalTime(), "Time")
             << " track.GetMomentumDirection() : " << track.GetMomentumDirection() << G4endl;
    }
#endif
  return &fParticleChange;
}

G4double G4ChemReboundTransportation::AlongStepGetPhysicalInteractionLength(
  const G4Track& track, G4double /*previousStepSize*/, G4double /*currentMinimumStep*/,
  G4double& /*currentSafety*/, G4GPILSelection* /*selection*/)
{
  if (!fpBoundingBox->contains(track.GetPosition())) {
    G4ExceptionDescription errMsg;
    errMsg << "Point is out of box : " << *fpBoundingBox
           << " of particle : " << GetIT(track)->GetName() << "(" << track.GetTrackID()
           << ") : " << track.GetPosition();
    G4Exception(
      "ChemReboundTransportation::AlongStepGetPhysicalInteractionLength"
      "ChemReboundTransportation",
      "ChemReboundTransportation", FatalErrorInArgument, errMsg);
  }
  if (fNistWater != track.GetMaterial()) {
    G4ExceptionDescription errMsg;
    errMsg << "This is not water";
    G4Exception(
      "ChemReboundTransportation::AlongStepGetPhysicalInteractionLength"
      "ChemReboundTransportation",
      "ChemReboundTransportation", FatalErrorInArgument, errMsg);
  }

  G4double geometryStepLength = DBL_MAX;
  State(theInteractionTimeLeft) = DBL_MAX;

  auto molConf = GetMolecule(track)->GetMolecularConfiguration();
  G4ITReactionPerTime& reactionPerTime = fReactionSet->GetReactionsPerTime();
  auto reaction_i = reactionPerTime.begin();
  if (reaction_i == reactionPerTime.end()) {
    State(fGeometryLimitedStep) = false;
    State(theInteractionTimeLeft) = fMaximumTimeStep;
    if (fVerboseLevel > 1) {
      G4cout << "out of reaction " << G4BestUnit(State(theInteractionTimeLeft), "Time") << G4endl;
    }
  }
  else {
    G4Track* pTrackA = (*reaction_i)->GetReactants().first;
    G4Track* pTrackB = (*reaction_i)->GetReactant(pTrackA);

    if (&track == pTrackA || &track == pTrackB) {
      State(theInteractionTimeLeft) = GetTimeToBoundary(track);
      State(fTimeStepReachedLimit) = false;
      State(fGeometryLimitedStep) = false;

      if (fVerboseLevel > 1)
        G4cout << "Molecule A is of type : " << GetMolecule(track)->GetName()
               << " with trackID : " << track.GetTrackID()
               << "  fMaximumTimeStep : " << G4BestUnit(fMaximumTimeStep, "Time")
               << "  State(theInteractionTimeLeft) : "
               << G4BestUnit(State(theInteractionTimeLeft), "Time") << G4endl;

      if (State(theInteractionTimeLeft) < fInternalMinTimeStep) {
        State(fTimeStepReachedLimit) = true;
        State(theInteractionTimeLeft) = fInternalMinTimeStep;
      }
      else if (State(theInteractionTimeLeft) > fMaximumTimeStep) {
        State(fTimeStepReachedLimit) = true;
        State(theInteractionTimeLeft) = fMaximumTimeStep;
      }
    }
    else {
      State(fGeometryLimitedStep) = false;
      State(theInteractionTimeLeft) = DBL_MAX;
    }
  }

  geometryStepLength = calculateDistanceFromTimeStep(molConf, State(theInteractionTimeLeft));
  State(fTransportEndPosition) =
    geometryStepLength * track.GetMomentumDirection() + track.GetPosition();
  State(fTimeStepReachedLimit) = true;
  State(fCandidateEndGlobalTime) = track.GetGlobalTime() + State(theInteractionTimeLeft);
  State(fEndGlobalTimeComputed) = true;

#ifdef G4VERBOSE
  if (fVerboseLevel > 1) {
    G4cout << "ChemReboundTransportation::AlongStepGetPhysicalInteractionLength = "
           << G4BestUnit(geometryStepLength, "Length") << "  "
           << G4BestUnit(State(theInteractionTimeLeft), "Time")
           << " | trackID = " << track.GetTrackID() << G4endl;
  }
#endif
  return geometryStepLength;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4ChemReboundTransportation::AlongStepDoIt(const G4Track& track,
                                                              const G4Step& step)
{
  if (GetIT(track)->GetTrackingInfo()->IsLeadingStep()) {
    G4double spaceStep = DBL_MAX;
    auto molConf = GetMolecule(track)->GetMolecularConfiguration();
    spaceStep = calculateDistanceFromTimeStep(molConf, State(theInteractionTimeLeft));

    State(fGeometryLimitedStep) = false;
    State(fTransportEndPosition) =
      BouncingAction(track.GetPosition() + spaceStep * G4RandomDirection());
    State(fEndPointDistance) = spaceStep;
    if (fVerboseLevel > 1)
    // if(GetMolecule(track)->GetName() == "e_aq^-1")

    {
      G4cout << "ChemReboundTransportation::AlongStepDoIt() :" << " trackID : "
             << track.GetTrackID()
             << " Molecule name: " << "prePosition : " << step.GetPreStepPoint()->GetPosition()
             << "  postPostion : " << step.GetPostStepPoint()->GetPosition() << "  "
             << GetMolecule(track)->GetName() << " State(theInteractionTimeLeft)  : "
             << G4BestUnit(State(theInteractionTimeLeft), "Time")
             << "   Diffusion length : " << G4BestUnit(step.GetStepLength(), "Length")
             << " within time step : " << G4BestUnit(step.GetDeltaTime(), "Time")
             << "\t Current global time : " << G4BestUnit(track.GetGlobalTime(), "Time")
             << "  track.GetMomentumDirection() : " << track.GetMomentumDirection() << G4endl;
    }
  }

  G4ITTransportation::AlongStepDoIt(track, step);

  return &fParticleChange;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector G4ChemReboundTransportation::BouncingAction(const G4ThreeVector& nextPosition)
{
  // from Karamitros, Mathieu et al.2020,arXiv:2006.14225 (2020)
  // https://doi.org/10.48550/arXiv.2006.14225

  G4ThreeVector output;

  G4double RxM = fpBoundingBox->Getxhi();
  G4double RyM = fpBoundingBox->Getyhi();
  G4double RzM = fpBoundingBox->Getzhi();
  G4double Rxm = fpBoundingBox->Getxlo();
  G4double Rym = fpBoundingBox->Getylo();
  G4double Rzm = fpBoundingBox->Getzlo();

  G4double x = calculateNextCoordinate(nextPosition.getX(), RxM, Rxm);
  G4double y = calculateNextCoordinate(nextPosition.getY(), RyM, Rym);
  G4double z = calculateNextCoordinate(nextPosition.getZ(), RzM, Rzm);
  output.set(x, y, z);
  return output;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ChemReboundTransportation::calculateNextCoordinate(G4double nextPos, G4double high,
                                                              G4double low)
{
  // from Karamitros, Mathieu et al.2020,arXiv:2006.14225 (2020)

  G4double length = high - low;
  if (std::abs(length) < 1e-10) {
    return low;
  }

  G4double relativePos = std::abs(nextPos - low);
  if (!std::isfinite(relativePos)) {
    return low;  // no crash
  }

  G4double n = relativePos / length;
  if (!std::isfinite(n)) {
    return low;
  }

  G4double truncVal = std::floor(n);//n is already positive
  G4double h = truncVal;
  if(truncVal > 2.0){
    h = std::fmod(truncVal, 2.0);
  }

  G4double mod = relativePos;

  if(relativePos > length) {
    mod = std::fmod(relativePos, length);
  }

  return low + h * length + (1 - 2 * h) * std::abs(mod);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ChemReboundTransportation::calculateDistanceFromTimeStep(MolConf mol, G4double timeStep)
{
  G4double diffuCoeff = mol->GetDiffusionCoefficient();
  if (mol->GetDiffusionCoefficient() <= 0) {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "GetDiffusionCoefficient is negative";
    G4Exception("ChemReboundTransportation::calculateDistanceFromTimeStep",
                "ChemReboundTransportation030", FatalErrorInArgument, exceptionDescription);
  }
  G4double sqrt_2Dt = sqrt(2 * diffuCoeff * timeStep);
  G4double x = G4RandGauss::shoot(0, sqrt_2Dt);
  G4double y = G4RandGauss::shoot(0, sqrt_2Dt);
  G4double z = G4RandGauss::shoot(0, sqrt_2Dt);
  return sqrt(x * x + y * y + z * z);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ChemReboundTransportation::GetTimeToBoundary(const G4Track& track)
{
  if (!fpBoundingBox->contains(track.GetPosition())) {
    G4ExceptionDescription errMsg;
    errMsg << "Point is out of box : " << *fpBoundingBox
           << " of particle : " << GetIT(track)->GetName() << "(" << track.GetTrackID()
           << ") : " << track.GetPosition();
    G4Exception(
      "BoundedBrownianAction::GetTimeToBoundary"
      "BoundedBrownianAction",
      "BoundedBrownianAction", FatalErrorInArgument, errMsg);
  }
  auto diffusionCoefficient = GetMolecule(track)->GetDiffusionCoefficient();

  auto dx = std::min(track.GetPosition().getX() - fpBoundingBox->Getxlo(),
                     fpBoundingBox->Getxhi() - track.GetPosition().getX());
  auto dy = std::min(track.GetPosition().getY() - fpBoundingBox->Getylo(),
                     fpBoundingBox->Getyhi() - track.GetPosition().getY());
  auto dz = std::min(track.GetPosition().getZ() - fpBoundingBox->Getzlo(),
                     fpBoundingBox->Getzhi() - track.GetPosition().getZ());

  std::vector<G4double> distanceVector{dx, dy, dz};
  G4double MinTime = DBL_MAX;
  for (const auto& it : distanceVector) {
    G4double distance = it;
    auto random = G4UniformRand();
    auto minTime = 1 / (4 * diffusionCoefficient) * pow(distance / InvErfc(random), 2);

    if (MinTime > minTime) {
      MinTime = minTime;
    }
  }
  return MinTime;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

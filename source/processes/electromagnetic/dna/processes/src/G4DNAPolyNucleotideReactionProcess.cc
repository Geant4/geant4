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
/// file: G4DNAPolyNucleotideReactionProcess.cc
/// brief: This file handls reaction process with DNA geometry.
#include "G4DNAPolyNucleotideReactionProcess.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4Molecule.hh"
#include "G4ITReactionChange.hh"
#include "G4IRTUtils.hh"
#include "G4VDNAHitModel.hh"
#include "G4LowEnergyEmProcessSubType.hh"

#ifndef PrepareState
#  define PrepareState()                                                       \
    G4PolyNucleotideReactionState* _state =                                    \
      this->GetState<G4PolyNucleotideReactionState>()
#endif

#ifndef State
#  define State(theXInfo) (_state->theXInfo)
#endif

G4DNAPolyNucleotideReactionProcess::G4DNAPolyNucleotideReactionProcess(
  const G4String& aName, G4int verbosityLevel)
  : G4VITDiscreteProcess(aName, fUserDefined)
  , fHasAlreadyReachedNullTime(false)
  , fVerbose(verbosityLevel)
  , fRCutOff(G4IRTUtils::GetDNADistanceCutOff())
  , fpDamageModel(nullptr)
{
  pParticleChange     = &fParticleChange;
  enableAtRestDoIt    = false;
  enableAlongStepDoIt = false;
  enablePostStepDoIt  = true;
  fProposesTimeStep   = true;
  SetProcessSubType(fLowEnergyStaticMol);
  G4VITProcess::SetInstantiateProcessState(false);
}

G4DNAPolyNucleotideReactionProcess::~G4DNAPolyNucleotideReactionProcess()
{
  delete fpDamageModel;  // for now, this process handls only one model
}

G4DNAPolyNucleotideReactionProcess::G4PolyNucleotideReactionState::
G4PolyNucleotideReactionState()
{
  fSampledMinTimeStep         = 0;
  fPreviousTimeAtPreStepPoint = -1;
}

G4double G4DNAPolyNucleotideReactionProcess::CalculateTimeStep(
  const G4Track& trackA, const G4double& /*userTimeStep*/)
{
  PrepareState();
  fHasAlreadyReachedNullTime      = false;
  State(fSampledMinTimeStep)      = DBL_MAX;
  State(theInteractionTimeLeft)   = DBL_MAX;
  State(currentInteractionLength) = -1;
#ifdef G4VERBOSE
  if(fVerbose > 1)
  {
    auto pMoleculeA = GetMolecule(trackA);
    G4cout << "________________________________________________________________"
              "_______"
           << G4endl;
    G4cout << "G4DNAPolyNucleotideReactionProcess::CalculateTimleStep"
           << G4endl;
    G4cout << "Check done for molecule : " << pMoleculeA->GetName() << " ("
           << trackA.GetTrackID() << ") " << G4endl;
  }
#endif
  //__________________________________________________________________
  // Retrieve general informations for making reactions
  G4double reactionTime =
    fpDamageModel->CalculateReactionTime(trackA, State(fNodeReactant));

  if(reactionTime < 0)
  {
    return DBL_MAX;
  }
  State(fSampledMinTimeStep)      = reactionTime;
  State(theInteractionTimeLeft)   = State(fSampledMinTimeStep);
  State(currentInteractionLength) = State(theInteractionTimeLeft);

#ifdef G4VERBOSE
  if(fVerbose > 1)
  {
    G4cout << " theInteractionTimeLeft : " << State(theInteractionTimeLeft)
           << G4endl;
    G4cout << " State(fNodeReactant) : " << State(fNodeReactant).index()
           << G4endl;

    G4cout << "________________________________________________________________"
              "_______"
           << G4endl;
  }
#endif

  return State(fSampledMinTimeStep);
}

G4double
G4DNAPolyNucleotideReactionProcess::PostStepGetPhysicalInteractionLength(
  const G4Track& track,
  G4double,  // previousStepSize
  G4ForceCondition* pForceCond)
{
  PrepareState();
  CalculateTimeStep(track);
  *pForceCond = NotForced;

  G4double previousTimeStep(-1.);

  if(State(fPreviousTimeAtPreStepPoint) != -1)
  {
    previousTimeStep =
      track.GetGlobalTime() - State(fPreviousTimeAtPreStepPoint);
  }

  State(fPreviousTimeAtPreStepPoint) = track.GetGlobalTime();

  if((fpState->currentInteractionLength <= 0) || (previousTimeStep < 0.0) ||
     (fpState->theNumberOfInteractionLengthLeft <= 0.0))
  {
    ResetNumberOfInteractionLengthLeft();
  }
  else if(previousTimeStep > 0.0)
  {
    SubtractNumberOfInteractionLengthLeft(
      previousTimeStep);  // TODO need to check this
  }
  return State(theInteractionTimeLeft) * -1;
  // return value * -1;//should regard this !
}

///////////////////////////////////////////////////////////////////////////////

G4VParticleChange* G4DNAPolyNucleotideReactionProcess::PostStepDoIt(
  const G4Track& track, const G4Step&)
{
  PrepareState();
  auto isReacted = fpDamageModel->DoReaction(
    track, State(theInteractionTimeLeft), State(fNodeReactant));
  if(!isReacted)
  {
    // no reaction even this process is called
    fParticleChange.Initialize(track);
    return pParticleChange;
  }
  fParticleChange.Initialize(track);  // To initialise TouchableChange
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  State(fPreviousTimeAtPreStepPoint) = -1;
  return pParticleChange;
}
void G4DNAPolyNucleotideReactionProcess::StartTracking(G4Track* track)
{
  G4VProcess::StartTracking(track);
  G4VITProcess::fpState.reset(new G4PolyNucleotideReactionState());
  G4VITProcess::StartTracking(track);
}

#undef State
#undef PrepareState

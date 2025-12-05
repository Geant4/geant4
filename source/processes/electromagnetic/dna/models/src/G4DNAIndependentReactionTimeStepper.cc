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
// 20/2/2019
// Author: HoangTRAN

#include "G4DNAIndependentReactionTimeStepper.hh"

#include "G4ChemicalMoleculeFinder.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMakeReaction.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DiffusionControlledReactionModel.hh"
#include "G4IRTUtils.hh"
#include "G4ITReactionChange.hh"
#include "G4Molecule.hh"
#include "G4Scheduler.hh"
#include "G4UnitsTable.hh"
#include "G4VDNAReactionModel.hh"
#include "Randomize.hh"

#include <memory>
using namespace std;
using namespace CLHEP;

G4DNAIndependentReactionTimeStepper::Utils::Utils(const G4Track& trackA, const G4Track& trackB)
  : fpTrackA(const_cast<G4Track*>(&trackA)), fpTrackB(const_cast<G4Track*>(&trackB))
{
  fpMoleculeA = GetMolecule(trackA);
  fpMoleculeB = GetMolecule(trackB);
  fUserMinTimeStep = 1 * CLHEP::ps;
}

G4DNAIndependentReactionTimeStepper::G4DNAIndependentReactionTimeStepper()
{
  fReactionSet->SortByTime();
}

void G4DNAIndependentReactionTimeStepper::Prepare()
{
  //fVerbose = G4Scheduler::Instance()->GetVerbose();
  if (G4Scheduler::Instance()->IsInteractionStep()) {
    fReactionSet->CleanAllReaction();
    fIsInitialized = false;
    fSampledPositions.clear();
    fSecondaries.clear();
    InitializeForNewTrack();
  }
}

void G4DNAIndependentReactionTimeStepper::InitializeForNewTrack()
{
  fCheckedTracks.clear();
  BuildChemicalMoleculeFinder()
}

G4double G4DNAIndependentReactionTimeStepper::CalculateStep(const G4Track& trackA,
                                                            const G4double& /*userMinTimeStep*/)
{
  auto pMoleculeA = GetMolecule(trackA);
  fSampledMinTimeStep = DBL_MAX;
  fCheckedTracks.insert(trackA.GetTrackID());

#ifdef G4VERBOSE
  if (fVerbose > 1) {
    G4cout << "________________________________________________________________"
              "_______"
           << G4endl;
    G4cout << "G4DNAIndependentReactionTimeStepper::CalculateStep" << G4endl;
    G4cout << "Check done for molecule : " << pMoleculeA->GetName() << " (" << trackA.GetTrackID()
           << ") " << G4endl;
  }
#endif

  auto pMolConfA = pMoleculeA->GetMolecularConfiguration();

  const auto pReactantList = fMolecularReactionTable->CanReactWith(pMolConfA);

  if (pReactantList == nullptr) {
    if (fVerbose > 1) {
      G4ExceptionDescription msg;
      msg << "G4DNAIndependentReactionTimeStepper::CalculateStep will return infinity "
             "for the reaction because the molecule "
          << pMoleculeA->GetName() << " does not have any reactants given in the reaction table."
          << G4endl;
      G4Exception("G4DNAIndependentReactionTimeStepper::CalculateStep",
                  "G4DNAIndependentReactionTimeStepper03", JustWarning, msg);
    }
    return DBL_MAX;
  }

  auto nbReactives = (G4int)pReactantList->size();

  if (nbReactives == 0) {
    if (fVerbose > 1) {
      G4ExceptionDescription msg;
      msg << "G4DNAIndependentReactionTimeStepper::CalculateStep will "
             "return infinity "
             "for the reaction because the molecule "
          << pMoleculeA->GetName() << " does not have any reactants given in the reaction table."
          << "This message can also result from a wrong implementation of "
             "the reaction table."
          << G4endl;
      G4Exception("G4DNAIndependentReactionTimeStepper::CalculateStep",
                  "G4DNAIndependentReactionTimeStepper04", JustWarning, msg);
    }
    return DBL_MAX;
  }
  fReactionModel->Initialise(pMolConfA, trackA);
  for (G4int i = 0; i < nbReactives; ++i) {
    auto pMoleculeB = (*pReactantList)[i];
    G4int key = pMoleculeB->GetMoleculeID();
    fRCutOff = G4IRTUtils::GetRCutOff();
    //______________________________________________________________
    // Retrieve reaction range
    const G4double Reff = fReactionModel->GetReactionRadius(i);
    std::vector<std::pair<G4TrackList::iterator, G4double>> resultIndices;
    resultIndices.clear();
    G4ChemicalMoleculeFinder::Instance()->FindNearestInRange(trackA, key, fRCutOff, resultIndices);

    if (resultIndices.empty()) {
      continue;
    }
    for (auto& it : resultIndices) {
      G4Track* pTrackB = *(std::get<0>(it));

      if (pTrackB == &trackA) {
        continue;
      }
      if (pTrackB == nullptr) {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription << "No trackB no valid";
        G4Exception(
          "G4DNAIndependentReactionTimeStepper"
          "::CalculateStep()",
          "G4DNAIndependentReactionTimeStepper007", FatalException, exceptionDescription);
      }
      if (fCheckedTracks.find(pTrackB->GetTrackID()) != fCheckedTracks.end()) {
        continue;
      }

      Utils utils(trackA, *pTrackB);
      auto pMolB = GetMolecule(pTrackB);
      auto pMolConfB = pMolB->GetMolecularConfiguration();
      G4double distance = (trackA.GetPosition() - pTrackB->GetPosition()).mag();
      if (distance * distance < Reff * Reff) {
        auto reactionData = fMolecularReactionTable->GetReactionData(pMolConfA, pMolConfB);
        if (G4Scheduler::Instance()->GetGlobalTime() == G4Scheduler::Instance()->GetStartTime()) {
          if (reactionData->GetProbability() > G4UniformRand()) {
            fSampledMinTimeStep = 0.;
          }
        }
      }
      else {
        G4double tempMinET = GetTimeToEncounter(trackA, *pTrackB);
        if (tempMinET < 0 || tempMinET > G4Scheduler::Instance()->GetEndTime()) {
          continue;
        }
        fSampledMinTimeStep = tempMinET;
        // Fixed the IRT_syn model (Stepper) to ensure it does not use minTimeStep (default = 1 ps)
        // when DNA reactions do not yet share
        // the same minTimeStep(The MinTimeStep is used to optimize the chemistry).
        // TODO: full test

        //if (tempMinET < fUserMinTimeStep) {
        //  fSampledMinTimeStep = fUserMinTimeStep;
        //}
      }
      CheckAndRecordResults(fSampledMinTimeStep, utils);
    }
  }
  return fSampledMinTimeStep;
}

void G4DNAIndependentReactionTimeStepper::CheckAndRecordResults(G4double reactionTime,
                                                                const Utils& utils)
{
  if (utils.fpTrackB->GetTrackStatus() != fAlive) {
    return;
  }

  if (&utils.fpTrackB == &utils.fpTrackA) {
    G4ExceptionDescription msg;
    msg << "A track is reacting with itself"
           " (which is impossible) ie fpTrackA == trackB"
        << G4endl;
    msg << "Molecule A is of type : " << utils.fpMoleculeA->GetName()
        << " with trackID : " << utils.fpTrackA->GetTrackID()
        << " and B : " << utils.fpMoleculeB->GetName()
        << " with trackID : " << utils.fpTrackB->GetTrackID() << G4endl;
    G4Exception("G4DNAIndependentReactionTimeStepper::RetrieveResults",
                "G4DNAIndependentReactionTimeStepper003", FatalErrorInArgument, msg);
  }

  if (fabs(utils.fpTrackB->GetGlobalTime() - utils.fpTrackA->GetGlobalTime())
      > utils.fpTrackA->GetGlobalTime() * (1. - 1. / 100))
  {
    // DEBUG
    G4ExceptionDescription msg;
    msg << "The interacting tracks are not synchronized in time" << G4endl;
    msg << "trackB->GetGlobalTime() != fpTrackA.GetGlobalTime()" << G4endl;

    msg << "fpTrackA : trackID : " << utils.fpTrackA->GetTrackID()
        << "\t Name :" << utils.fpMoleculeA->GetName()
        << "\t fpTrackA->GetGlobalTime() = " << G4BestUnit(utils.fpTrackA->GetGlobalTime(), "Time")
        << G4endl;

    msg << "trackB : trackID : " << utils.fpTrackB->GetTrackID()
        << "\t Name :" << utils.fpMoleculeB->GetName()
        << "\t trackB->GetGlobalTime() = " << G4BestUnit(utils.fpTrackB->GetGlobalTime(), "Time")
        << G4endl;

    G4Exception("G4DNAIndependentReactionTimeStepper::RetrieveResults",
                "G4DNAIndependentReactionTimeStepper004", FatalErrorInArgument, msg);
  }
  if (reactionTime < 0) {
    // DEBUG
    G4ExceptionDescription msg;
    msg << "The interacting tracks are not in good time" << G4endl;

    msg << "fpTrackA : trackID : " << utils.fpTrackA->GetTrackID()
        << "\t Name :" << utils.fpMoleculeA->GetName()
        << "\t fpTrackA->GetGlobalTime() = " << G4BestUnit(utils.fpTrackA->GetGlobalTime(), "Time")
        << G4endl;

    msg << "trackB : trackID : " << utils.fpTrackB->GetTrackID()
        << "\t Name :" << utils.fpMoleculeB->GetName()
        << "\t trackB->GetGlobalTime() = " << G4BestUnit(utils.fpTrackB->GetGlobalTime(), "Time")
        << G4endl;

    G4Exception("G4DNAIndependentReactionTimeStepper::CheckAndRecordResults",
                "G4DNAIndependentReactionTimeStepper1", FatalErrorInArgument, msg);
  }

  G4double globalTime = G4Scheduler::Instance()->GetGlobalTime();

  fReactionSet->AddReaction(reactionTime + globalTime, utils.fpTrackA, utils.fpTrackB);
  fSampledPositions[utils.fpTrackA->GetTrackID()] = utils.fpTrackA->GetPosition();
  fSampledPositions[utils.fpTrackB->GetTrackID()] = utils.fpTrackB->GetPosition();
}

std::unique_ptr<G4ITReactionChange> G4DNAIndependentReactionTimeStepper::FindReaction(
  G4ITReactionSet* pReactionSet, G4double& currentStepTime, const G4double globalTime)
{
  if (pReactionSet == nullptr) {
    return nullptr;
  }

  G4ITReactionPerTime& reactionPerTime = pReactionSet->GetReactionsPerTime();
  if (reactionPerTime.empty()) {
    return nullptr;
  }
  for (auto reaction_i = reactionPerTime.begin(); reaction_i != reactionPerTime.end();
       reaction_i = reactionPerTime.begin())
  {
    G4Track* pTrackA = (*reaction_i)->GetReactants().first;
    currentStepTime = DBL_MAX;
    if (pTrackA->GetTrackStatus() == fStopAndKill) {
      continue;
    }
    G4Track* pTrackB = (*reaction_i)->GetReactant(pTrackA);
    currentStepTime = DBL_MAX;
    if (pTrackB->GetTrackStatus() == fStopAndKill) {
      continue;
    }

    if (pTrackB == pTrackA) {
      G4ExceptionDescription msg;
      msg << "The IT reaction process sent back a reaction "
             "between trackA and trackB. ";
      msg << "The problem is trackA == trackB";
      G4Exception("G4DNAIndependentReactionTimeStepper::FindReaction",
                  "G4DNAIndependentReactionTimeStepper02", FatalErrorInArgument, msg);
    }
    G4double reactionTime = (*reaction_i)->GetTime();
    currentStepTime =  reactionTime - globalTime;
    if(fVerbose > 1)
    G4cout << " reaction Time : " << reactionTime << "  currentStepTime : " << currentStepTime
           << " globalTime : " << globalTime << "  " << pTrackA->GetTrackID() << " + "
           << pTrackB->GetTrackID() << G4endl;

    pReactionSet->SelectThisReaction(*reaction_i);
    if (fpReactionProcess != nullptr) {
      if ((fSampledPositions.find(pTrackA->GetTrackID()) == fSampledPositions.end()
           && (fSampledPositions.find(pTrackB->GetTrackID()) == fSampledPositions.end())))
      {
        G4ExceptionDescription msg;
        msg << "The positions of trackA and trackB have no counted ";
        G4Exception("G4DNAIndependentReactionTimeStepper::FindReaction",
                    "G4DNAIndependentReactionTimeStepper0001", FatalErrorInArgument, msg);
      }

      pTrackA->SetPosition(fSampledPositions[pTrackA->GetTrackID()]);
      pTrackB->SetPosition(fSampledPositions[pTrackB->GetTrackID()]);
      auto pReactionChange = fpReactionProcess->MakeReaction(*pTrackA, *pTrackB);
      if (pReactionChange == nullptr) {
        return nullptr;
      }
      G4int nbSecondaries = pReactionChange->GetNumberOfSecondaries();
      if (nbSecondaries > 0) {
        const std::vector<G4Track*>* productsVector = pReactionChange->GetfSecondary();
        for (const auto& it : *productsVector) {
          fSecondaries.push_back(it);
        }
      }
      return pReactionChange;
    }
  }
  return nullptr;
}

void G4DNAIndependentReactionTimeStepper::SetReactionModel(G4VDNAReactionModel* pReactionModel)
{
  fReactionModel = pReactionModel;
}

G4VDNAReactionModel* G4DNAIndependentReactionTimeStepper::GetReactionModel()
{
  return fReactionModel;
}

void G4DNAIndependentReactionTimeStepper::SetVerbose(G4int flag)
{
  fVerbose = flag;
}

G4double G4DNAIndependentReactionTimeStepper::GetTimeToEncounter(const G4Track& trackA,
                                                                 const G4Track& trackB)
{
  G4double timeToReaction = dynamic_cast<G4DiffusionControlledReactionModel*>(fReactionModel)
                              ->GetTimeToEncounter(trackA, trackB);
  return timeToReaction;
}

void G4DNAIndependentReactionTimeStepper::SetReactionProcess(G4VITReactionProcess* pReactionProcess)
{
  fpReactionProcess = pReactionProcess;
}
G4double G4DNAIndependentReactionTimeStepper::CalculateMinTimeStep(G4double currentGlobalTime,
                                                                   G4double /*definedMinTimeStep*/)
{
  G4double fTSTimeStep = DBL_MAX;
  // fUserMinTimeStep = definedMinTimeStep;
  if (!fIsInitialized) {
    InitializeReactions(currentGlobalTime);
  }

  G4int nbPreviousSecondaries = (G4int)fSecondaries.size();
  if (nbPreviousSecondaries > 0) {
    InitializeForNewTrack();
    for (const auto& it : fSecondaries) {
      CalculateStep(*it, fUserMinTimeStep);
    }
    fSecondaries.clear();
  }
  fTSTimeStep = GetNextReactionTime() - currentGlobalTime;
  if (fTSTimeStep < 0) {
    G4ExceptionDescription msg;
    msg << "fTSTimeStep < 0" << ": fTSTimeStep : " << fTSTimeStep
        << " GetNextReactionTime() : " << GetNextReactionTime()
        << "  currentGlobalTime : " << currentGlobalTime << G4endl;
    G4Exception("G4DNAIndependentReactionTimeStepper::CalculateMinTimeStep",
                "G4DNAIndependentReactionTimeStepper002", FatalErrorInArgument, msg);
  }
  return fTSTimeStep;
}

void G4DNAIndependentReactionTimeStepper::InitializeReactions(G4double /*currentGlobalTime*/)
{
  fCheckedTracks.clear();
  for (auto pTrack : *fpTrackContainer->GetMainList()) {
    if (pTrack == nullptr) {
      G4ExceptionDescription msg;
      msg << "No track found.";
      G4Exception("G4DNAIndependentReactionTimeStepper::InitializeReactions",
                  "G4DNAIndependentReactionTimeStepper030", FatalErrorInArgument, msg);
      continue;
    }

    G4TrackStatus trackStatus = pTrack->GetTrackStatus();
    if (trackStatus == fStopAndKill || trackStatus == fStopButAlive) {
      continue;
    }
    CalculateStep(*pTrack, fUserMinTimeStep);
  }
  if (fVerbose > 0)
    G4cout << "InitializeReactions : reaction events : "
           << fReactionSet->GetReactionsPerTime().size() << ". The previous time step : "
           << G4BestUnit(G4Scheduler::Instance()->GetPreviousTimeStep(), "Time") << G4endl;
  fIsInitialized = true;
}

G4double G4DNAIndependentReactionTimeStepper::GetNextReactionTime()
{
  G4double output = DBL_MAX;
  auto nextReaction = GetNextReaction();
  if (nextReaction != nullptr) {
    output = GetNextReaction()->GetTime();
  }
  return output;
}

const G4ITReaction* G4DNAIndependentReactionTimeStepper::GetNextReaction()
{
  G4ITReaction* output = nullptr;
  G4ITReactionPerTime& reactionPerTime = fReactionSet->GetReactionsPerTime();
  auto reaction_i = reactionPerTime.begin();
  if (reaction_i != reactionPerTime.end()) {
    output = (reaction_i->get());
  }
  return output;
}
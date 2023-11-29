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

#include <memory>
#include "G4DNAIndependentReactionTimeStepper.hh"
#include "G4VDNAReactionModel.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4UnitsTable.hh"
#include "G4Molecule.hh"
#include "G4ChemicalMoleculeFinder.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMakeReaction.hh"
#include "G4ITReactionChange.hh"
#include "G4Scheduler.hh"
#include "G4IRTUtils.hh"
#include "Randomize.hh"
#include "G4DiffusionControlledReactionModel.hh"
using namespace std;
using namespace CLHEP;

G4DNAIndependentReactionTimeStepper::Utils::Utils(const G4Track& trackA,
                                                  const G4Track& trackB)
  : fTrackA(trackA)
  , fTrackB(trackB)
{
  fpMoleculeA = GetMolecule(trackA);
  fpMoleculeB = GetMolecule(trackA);
}

G4DNAIndependentReactionTimeStepper::G4DNAIndependentReactionTimeStepper()
  : G4VITTimeStepComputer()
{
  fReactionSet->SortByTime();
}

void G4DNAIndependentReactionTimeStepper::Prepare()
{
  G4VITTimeStepComputer::Prepare();
  fSampledPositions.clear();
  BuildChemicalMoleculeFinder()
}

void G4DNAIndependentReactionTimeStepper::InitializeForNewTrack()
{
  if(fReactants != nullptr)
  {
    fReactants.reset();
  }
  fSampledMinTimeStep        = DBL_MAX;
  fHasAlreadyReachedNullTime = false;
}

G4double G4DNAIndependentReactionTimeStepper::CalculateStep(
  const G4Track& trackA, const G4double& userMinTimeStep)
{
  auto pMoleculeA = GetMolecule(trackA);
  InitializeForNewTrack();
  fUserMinTimeStep = userMinTimeStep;
  fCheckedTracks.insert(trackA.GetTrackID());

#ifdef G4VERBOSE
  if(fVerbose != 0)
  {
    G4cout << "________________________________________________________________"
              "_______"
           << G4endl;
    G4cout << "G4DNAIndependentReactionTimeStepper::CalculateStep" << G4endl;
    G4cout << "Check done for molecule : " << pMoleculeA->GetName() << " ("
           << trackA.GetTrackID() << ") " << G4endl;
  }
#endif

  auto pMolConfA = pMoleculeA->GetMolecularConfiguration();

  const auto pReactantList = fMolecularReactionTable->CanReactWith(pMolConfA);

  if(pReactantList == nullptr)
  {
#ifdef G4VERBOSE
    if(fVerbose > 1)
    {
      G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
      G4cout << "!!! WARNING" << G4endl;
      G4cout << "G4DNAIndependentReactionTimeStepper::CalculateStep will "
                "return infinity "
                "for the reaction because the molecule "
             << pMoleculeA->GetName()
             << " does not have any reactants given in the reaction table."
             << G4endl;
      G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    }
#endif
    return DBL_MAX;
  }

  G4int nbReactives = (G4int)pReactantList->size();

  if(nbReactives == 0)
  {
#ifdef G4VERBOSE
    //    DEBUG
    if(fVerbose != 0)
    {
      G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
      G4cout << "!!! WARNING" << G4endl;
      G4cout << "G4DNAIndependentReactionTimeStepper::CalculateStep will "
                "return infinity "
                "for the reaction because the molecule "
             << pMoleculeA->GetName()
             << " does not have any reactants given in the reaction table."
             << "This message can also result from a wrong implementation of "
                "the reaction table."
             << G4endl;
      G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    }
#endif
    return DBL_MAX;
  }
  fReactants = std::make_shared<vector<G4Track*>>();
  fReactionModel->Initialise(pMolConfA, trackA);
  for(G4int i = 0; i < nbReactives; ++i)
  {
    auto pMoleculeB = (*pReactantList)[i];
    G4int key       = pMoleculeB->GetMoleculeID();

    // fRCutOff = G4IRTUtils::GetRCutOff(1 * ps);
    fRCutOff = G4IRTUtils::GetRCutOff();
    //______________________________________________________________
    // Retrieve reaction range
    const G4double Reff = fReactionModel->GetReactionRadius(i);
    std::vector<std::pair<G4TrackList::iterator, G4double>> resultIndices;
    resultIndices.clear();
    G4ChemicalMoleculeFinder::Instance()->FindNearestInRange(
      trackA, key, fRCutOff, resultIndices);

    if(resultIndices.empty())
    {
      continue;
    }
    for(auto& it : resultIndices)
    {
      G4Track* pTrackB = *(std::get<0>(it));

      if(pTrackB == &trackA)
      {
        continue;
      }
      if(pTrackB == nullptr)
      {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription << "No trackB no valid";
        G4Exception("G4DNAIndependentReactionTimeStepper"
                    "::CalculateStep()",
                    "G4DNAIndependentReactionTimeStepper007", FatalException,
                    exceptionDescription);
      }else
      {
        if(fCheckedTracks.find(pTrackB->GetTrackID()) != fCheckedTracks.end())
        {
          continue;
        }

        Utils utils(trackA, *pTrackB);

        auto pMolB        = GetMolecule(pTrackB);
        auto pMolConfB    = pMolB->GetMolecularConfiguration();
        G4double distance = (trackA.GetPosition() - pTrackB->GetPosition()).mag();
        if(distance * distance < Reff * Reff)
        {
          auto reactionData =
            fMolecularReactionTable->GetReactionData(pMolConfA, pMolConfB);
          if(G4Scheduler::Instance()->GetGlobalTime() == G4Scheduler::Instance()->GetStartTime())
          {
            if(reactionData->GetProbability() > G4UniformRand())
            {
              if(!fHasAlreadyReachedNullTime)
              {
                fReactants->clear();
                fHasAlreadyReachedNullTime = true;
              }
              fSampledMinTimeStep = 0.;
              CheckAndRecordResults(utils);
            }
          }
        }
        else
        {
          G4double tempMinET = GetTimeToEncounter(trackA, *pTrackB);
          if(tempMinET < 0 || tempMinET > G4Scheduler::Instance()->GetEndTime())
          {
            continue;
          }
          if(tempMinET >= fSampledMinTimeStep)
          {
            continue;
          }
          fSampledMinTimeStep = tempMinET;
          fReactants->clear();
          CheckAndRecordResults(utils);
        }
      }
    }
  }

#ifdef G4VERBOSE
  if(fVerbose != 0)
  {
    G4cout << "G4DNAIndependentReactionTimeStepper::CalculateStep will finally "
              "return :"
           << G4BestUnit(fSampledMinTimeStep, "Time") << G4endl;

    if(fVerbose > 1)
    {
      G4cout << "Selected reactants for trackA: " << pMoleculeA->GetName()
             << " (" << trackA.GetTrackID() << ") are: ";

      vector<G4Track*>::iterator it;
      for(it = fReactants->begin(); it != fReactants->end(); it++)
      {
        G4Track* trackB = *it;
        G4cout << GetMolecule(trackB)->GetName() << " (" << trackB->GetTrackID()
               << ") \t ";
      }
      G4cout << G4endl;
    }
  }
#endif
  return fSampledMinTimeStep;
}

void G4DNAIndependentReactionTimeStepper::CheckAndRecordResults(
  const Utils& utils)
{
  if(utils.fTrackB.GetTrackStatus() != fAlive)
  {
    return;
  }

  if(&utils.fTrackB == &utils.fTrackA)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "A track is reacting with itself"
                            " (which is impossible) ie fpTrackA == trackB"
                         << G4endl;
    exceptionDescription << "Molecule A is of type : "
                         << utils.fpMoleculeA->GetName()
                         << " with trackID : " << utils.fTrackA.GetTrackID()
                         << " and B : " << utils.fpMoleculeB->GetName()
                         << " with trackID : " << utils.fTrackB.GetTrackID()
                         << G4endl;
    G4Exception("G4DNAIndependentReactionTimeStepper::RetrieveResults",
                "G4DNAIndependentReactionTimeStepper003", FatalErrorInArgument,
                exceptionDescription);
  }

  if(fabs(utils.fTrackB.GetGlobalTime() - utils.fTrackA.GetGlobalTime()) >
     utils.fTrackA.GetGlobalTime() * (1. - 1. / 100))
  {
    // DEBUG
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
      << "The interacting tracks are not synchronized in time" << G4endl;
    exceptionDescription
      << "trackB->GetGlobalTime() != fpTrackA.GetGlobalTime()" << G4endl;

    exceptionDescription << "fpTrackA : trackID : "
                         << utils.fTrackA.GetTrackID()
                         << "\t Name :" << utils.fpMoleculeA->GetName()
                         << "\t fpTrackA->GetGlobalTime() = "
                         << G4BestUnit(utils.fTrackA.GetGlobalTime(), "Time")
                         << G4endl;

    exceptionDescription << "trackB : trackID : " << utils.fTrackB.GetTrackID()
                         << "\t Name :" << utils.fpMoleculeB->GetName()
                         << "\t trackB->GetGlobalTime() = "
                         << G4BestUnit(utils.fTrackB.GetGlobalTime(), "Time")
                         << G4endl;

    G4Exception("G4DNAIndependentReactionTimeStepper::RetrieveResults",
                "G4DNAIndependentReactionTimeStepper004", FatalErrorInArgument,
                exceptionDescription);
  }
  fReactants->push_back(const_cast<G4Track*>(&utils.fTrackB));
}

std::unique_ptr<G4ITReactionChange>
G4DNAIndependentReactionTimeStepper::FindReaction(
  G4ITReactionSet* pReactionSet, const G4double& currentStepTime,
  const G4double& /*previousStepTime*/,
  const G4bool& /*reachedUserStepTimeLimit*/)
{
  if(pReactionSet == nullptr)
  {
    return nullptr;
  }

  G4ITReactionPerTime& reactionPerTime = pReactionSet->GetReactionsPerTime();
  if(reactionPerTime.empty())
  {
    return nullptr;
  }

  for(auto reaction_i                                 = reactionPerTime.begin();
      reaction_i != reactionPerTime.end(); reaction_i = reactionPerTime.begin())
  {
    G4Track* pTrackA = (*reaction_i)->GetReactants().first;
    if(pTrackA->GetTrackStatus() == fStopAndKill)
    {
      continue;
    }
    G4Track* pTrackB = (*reaction_i)->GetReactant(pTrackA);
    if(pTrackB->GetTrackStatus() == fStopAndKill)
    {
      continue;
    }

    if(pTrackB == pTrackA)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "The IT reaction process sent back a reaction "
                              "between trackA and trackB. ";
      exceptionDescription << "The problem is trackA == trackB";
      G4Exception("G4DNAIndependentReactionTimeStepper::FindReaction",
                  "G4DNAIndependentReactionTimeStepper02", FatalErrorInArgument,
                  exceptionDescription);
    }
    pReactionSet->SelectThisReaction(*reaction_i);
    if(fpReactionProcess != nullptr &&
       fpReactionProcess->TestReactibility(*pTrackA, *pTrackB, currentStepTime,
                                           false))
    {
      if((fSampledPositions.find(pTrackA->GetTrackID()) ==
            fSampledPositions.end() &&
          (fSampledPositions.find(pTrackB->GetTrackID()) ==
           fSampledPositions.end())))
      {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription
          << "The positions of trackA and trackB have no counted ";
        G4Exception("G4DNAIndependentReactionTimeStepper::FindReaction",
                    "G4DNAIndependentReactionTimeStepper0001",
                    FatalErrorInArgument, exceptionDescription);
      }

      pTrackA->SetPosition(fSampledPositions[pTrackA->GetTrackID()]);
      pTrackB->SetPosition(fSampledPositions[pTrackB->GetTrackID()]);
      auto pReactionChange =
        fpReactionProcess->MakeReaction(*pTrackA, *pTrackB);
      if(pReactionChange == nullptr)
      {
        return nullptr;
      }
      return pReactionChange;
    }
  }
  return nullptr;
}

void G4DNAIndependentReactionTimeStepper::SetReactionModel(
  G4VDNAReactionModel* pReactionModel)
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

G4double G4DNAIndependentReactionTimeStepper::GetTimeToEncounter(
  const G4Track& trackA, const G4Track& trackB)
{
  G4double timeToReaction =
    dynamic_cast<G4DiffusionControlledReactionModel*>(fReactionModel)
      ->GetTimeToEncounter(trackA, trackB);
  return timeToReaction;
}

void G4DNAIndependentReactionTimeStepper::SetReactionProcess(
  G4VITReactionProcess* pReactionProcess)
{
  fpReactionProcess = pReactionProcess;
}
G4double G4DNAIndependentReactionTimeStepper::CalculateMinTimeStep(
  G4double /*currentGlobalTime*/, G4double definedMinTimeStep)
{
  G4double fTSTimeStep = DBL_MAX;
  fCheckedTracks.clear();

  for(auto pTrack : *fpTrackContainer->GetMainList())
  {
    if(pTrack == nullptr)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "No track found.";
      G4Exception("G4DNAIndependentReactionTimeStepper::CalculateMinTimeStep",
                  "G4DNAIndependentReactionTimeStepper006",
                  FatalErrorInArgument, exceptionDescription);
      continue;
    }

    G4TrackStatus trackStatus = pTrack->GetTrackStatus();
    if(trackStatus == fStopAndKill || trackStatus == fStopButAlive)
    {
      continue;
    }

    G4double sampledMinTimeStep   = CalculateStep(*pTrack, definedMinTimeStep);
    G4TrackVectorHandle reactants = GetReactants();

    if(sampledMinTimeStep < fTSTimeStep)
    {
      fTSTimeStep = sampledMinTimeStep;
      fReactionSet->CleanAllReaction();
      if(reactants)
      {
        fReactionSet->AddReactions(fTSTimeStep, const_cast<G4Track*>(pTrack),
                                   reactants);

        fSampledPositions[pTrack->GetTrackID()] = pTrack->GetPosition();
        for(const auto& it : *fReactants)
        {
          auto pTrackB = it;
          // G4cout<<"position : "<<pTrackB->GetTrackID()<<G4endl;
          fSampledPositions[pTrackB->GetTrackID()] = pTrackB->GetPosition();
        }
        ResetReactants();
      }
    }
    else if(fTSTimeStep == sampledMinTimeStep && G4bool(reactants))
    {
      fReactionSet->AddReactions(fTSTimeStep, const_cast<G4Track*>(pTrack),
                                 reactants);

      fSampledPositions[pTrack->GetTrackID()] = pTrack->GetPosition();
      for(const auto& it : *fReactants)
      {
        auto pTrackB = it;
        // G4cout<<"position : "<<pTrackB->GetTrackID()<<G4endl;
        fSampledPositions[pTrackB->GetTrackID()] = pTrackB->GetPosition();
      }
      ResetReactants();
    }
    else if(reactants)
    {
      ResetReactants();
    }
  }

  return fTSTimeStep;
}

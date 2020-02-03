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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITModelProcessor.hh"
#include "G4VITTimeStepComputer.hh"
#include "G4VITReactionProcess.hh"
#include "G4ITReaction.hh"
#include "G4ITTrackHolder.hh"
#include "G4ITTrackingManager.hh"
#include "G4VITStepModel.hh"
#include "G4UserTimeStepAction.hh"
#include "G4UnitsTable.hh"
#include <vector>

//#define DEBUG_MEM

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
#endif

G4ITModelProcessor::G4ITModelProcessor()
{
    fpTrack = nullptr;
    fInitialized = false;
    fUserMinTimeStep = -1.;
    fTSTimeStep = DBL_MAX;
    fpTrackingManager = nullptr;
    fReactionSet = nullptr;
    fpTrackContainer = nullptr;
    fpModelHandler = nullptr;
    fpActiveModelWithMinTimeStep = nullptr;
    fComputeTimeStep = false;
    fComputeReaction = false;
}

G4ITModelProcessor::~G4ITModelProcessor() = default;

void G4ITModelProcessor::RegisterModel(double time, G4VITStepModel* model)
{
    fpModelHandler->RegisterModel(model, time);
}

void G4ITModelProcessor::Initialize()
{
    fpModelHandler->Initialize();
    fReactionSet = G4ITReactionSet::Instance();
    fpTrackContainer = G4ITTrackHolder::Instance();
    fInitialized = true;
    fComputeTimeStep = false;
    fComputeReaction = false;
    if (fpModelHandler->GetTimeStepComputerFlag())
    {
        fComputeTimeStep = true;
    }
    if (fpModelHandler->GetReactionProcessFlag())
    {
        fComputeReaction = true;
    }
}

G4double G4ITModelProcessor::CalculateMinTimeStep(G4double currentGlobalTime,
                                                  G4double definedMinTimeStep)
{

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
    MemStat mem_first, mem_second, mem_diff;
    mem_first = MemoryUsage();
#endif

    fpActiveModelWithMinTimeStep = nullptr;
    fTSTimeStep = DBL_MAX;

    InitializeStepper(currentGlobalTime, definedMinTimeStep);

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
    mem_second = MemoryUsage();
    mem_diff = mem_second-mem_first;
    G4cout << "\t || MEM || G4Scheduler::CalculateMinTimeStep || After "
    "computing fpModelProcessor -> InitializeStepper, diff is : "
    << mem_diff
    << G4endl;
#endif

    for (auto pTrack : *fpTrackContainer->GetMainList())
    {
        if (pTrack == nullptr)
        {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription << "No track found.";
            G4Exception("G4Scheduler::CalculateMinStep", "ITScheduler006",
                        FatalErrorInArgument, exceptionDescription);
            continue;
        }

#ifdef DEBUG
        G4cout << "*_* " << GetIT(track)->GetName()
        << " ID: " << track->GetTrackID()
        << " at time : " << track->GetGlobalTime()
        << G4endl;
#endif

        G4TrackStatus trackStatus = pTrack->GetTrackStatus();
        if (trackStatus == fStopAndKill || trackStatus == fStopButAlive)
        {
            continue;
        }

        CalculateTimeStep(pTrack, definedMinTimeStep);
        // if MT mode at track level, this command should be displaced
        ExtractTimeStepperData();
    }

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
    mem_second = MemoryUsage();
    mem_diff = mem_second-mem_first;
    G4cout << "\t || MEM || G4Scheduler::CalculateMinTimeStep || "
    "After looping on tracks, diff is : " << mem_diff << G4endl;
#endif
    return fTSTimeStep;
}

//_________________________________________________________________________

void G4ITModelProcessor::ExtractTimeStepperData()
{
    if (fpTrack == nullptr)
    {
        CleanProcessor();
        return;
    }

    for (auto pStepModel : fActiveModels)
    {
        if (pStepModel == nullptr)
        {
            continue;
        }

        auto pTimeStepper = pStepModel->GetTimeStepper();
        G4double sampledMinTimeStep = pTimeStepper->GetSampledMinTimeStep();
        G4TrackVectorHandle reactants = pTimeStepper->GetReactants();

        if (sampledMinTimeStep < fTSTimeStep)
        {
            fpActiveModelWithMinTimeStep = pStepModel;
            fTSTimeStep = sampledMinTimeStep;
            //fReactingTracks.clear();

            fReactionSet->CleanAllReaction();
            if (reactants)
            {
                //  fReactingTracks.insert(make_pair(track, reactants));
                fReactionSet->AddReactions(fTSTimeStep,
                                           const_cast<G4Track*>(fpTrack),
                                           reactants);
                pTimeStepper->ResetReactants();
            }
        }
        else if (fTSTimeStep == sampledMinTimeStep && bool(reactants))
        {
            // fReactingTracks.insert(make_pair(track, reactants));
            fReactionSet->AddReactions(fTSTimeStep,
                                       const_cast<G4Track*>(fpTrack),
                                       reactants);
            pTimeStepper->ResetReactants();
        }
        else if (reactants)
        {
            pTimeStepper->ResetReactants();
        }
    }

    CleanProcessor();
}

//______________________________________________________________________________

void G4ITModelProcessor::InitializeStepper(G4double currentGlobalTime,
                                           G4double userMinTime)
{
    G4VITTimeStepComputer::SetTimes(currentGlobalTime, userMinTime);

#if defined (DEBUG_MEM)
    MemStat mem_first, mem_second, mem_diff;
            mem_first = MemoryUsage();
#endif

    fActiveModels = fpModelHandler->GetActiveModels(currentGlobalTime);

    for (auto& pModel : fActiveModels)
    {
        pModel->PrepareNewTimeStep();
    }

#if defined (DEBUG_MEM)
    mem_second = MemoryUsage();
            mem_diff = mem_second-mem_first;
            G4cout << "\t || MEM || G4ITModelProcessor::InitializeStepper || After computing stepper -> Prepare(), diff is : " << mem_diff << G4endl;
#endif

}

//______________________________________________________________________________
void G4ITModelProcessor::CalculateTimeStep(const G4Track* pTrack,
                                           const G4double userMinTimeStep)
{
    CleanProcessor();
    if (pTrack == nullptr)
    {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription << "No track was passed to the method.";
        G4Exception("G4ITModelProcessor::CalculateStep",
                    "ITModelProcessor004",
                    FatalErrorInArgument,
                    exceptionDescription);
    }
    SetTrack(pTrack);
    fUserMinTimeStep = userMinTimeStep;

    DoCalculateStep();
}

//______________________________________________________________________________

void G4ITModelProcessor::DoCalculateStep()
{
    for (auto& pStepModel : fActiveModels)
    {
        pStepModel->GetTimeStepper()->CalculateStep(*fpTrack, fUserMinTimeStep);
    }
}

//_________________________________________________________________________

void G4ITModelProcessor::ComputeTrackReaction(G4ITStepStatus fITStepStatus,
                                              G4double fGlobalTime,
                                              G4double currentTimeStep,
                                              G4double previousTimeStep,
                                              G4bool reachedUserTimeLimit,
                                              G4double fTimeTolerance,
                                              G4UserTimeStepAction* fpUserTimeStepAction,
                                              G4int
#ifdef G4VERBOSE
fVerbose
#endif
)
{
//  if (fReactingTracks.empty())
    if (fReactionSet->Empty())
    {
        return;
    }

    if (fITStepStatus == eCollisionBetweenTracks)
        //        if(fInteractionStep == false)
    {
        // TODO
        FindReaction(fReactionSet,
                     currentTimeStep,
                     previousTimeStep,
                     reachedUserTimeLimit);
        // TODO
        // A ne faire uniquement si le temps choisis est celui calculé par le time stepper
        // Sinon utiliser quelque chose comme : fModelProcessor->FindReaction(&fMainList);

        for (auto& pChanges : fReactionInfo)
        {
            auto pTrackA = const_cast<G4Track*>(pChanges->GetTrackA());
            auto pTrackB = const_cast<G4Track*>(pChanges->GetTrackB());

            if (pTrackA == nullptr
                || pTrackB == nullptr
                || pTrackA->GetTrackStatus() == fStopAndKill
                || pTrackB->GetTrackStatus() == fStopAndKill)
            {
                continue;
            }

            G4int nbSecondaries = pChanges->GetNumberOfSecondaries();
            const std::vector<G4Track*>* productsVector = pChanges->GetfSecondary();

            if (fpUserTimeStepAction)
            {
                fpUserTimeStepAction->UserReactionAction(*pTrackA,
                                                         *pTrackB,
                                                         productsVector);
            }

#ifdef G4VERBOSE
            if (fVerbose)
            {
                G4cout << "At time : " << std::setw(7) << G4BestUnit(fGlobalTime, "Time")
                       << " Reaction : " << GetIT(pTrackA)->GetName() << " ("
                       << pTrackA->GetTrackID() << ") + " << GetIT(pTrackB)->GetName() << " ("
                       << pTrackB->GetTrackID() << ") -> ";
            }
#endif

            if (nbSecondaries > 0)
            {
                for (int i = 0; i < nbSecondaries; ++i)
                {
#ifdef G4VERBOSE
                    if (fVerbose && i != 0)
                    {
                        G4cout << " + ";
                    }
#endif

                    G4Track* secondary = (*productsVector)[i]; //changes->GetSecondary(i);
                    fpTrackContainer->_PushTrack(secondary);
                    GetIT(secondary)->SetParentID(pTrackA->GetTrackID(),
                                                  pTrackB->GetTrackID());

                    if (secondary->GetGlobalTime() - fGlobalTime > fTimeTolerance)
                    {
                        G4ExceptionDescription exceptionDescription;
                        exceptionDescription << "The time of the secondary should not be bigger than the"
                                                " current global time."
                                             << " This may cause synchronization problem. If the process you"
                                                " are using required "
                                             << "such feature please contact the developpers." << G4endl
                                             << "The global time in the step manager : "
                                             << G4BestUnit(fGlobalTime, "Time")
                                             << G4endl
                                             << "The global time of the track : "
                                             << G4BestUnit(secondary->GetGlobalTime(), "Time")
                                             << G4endl;

                        G4Exception("G4Scheduler::ComputeInteractionBetweenTracks",
                                    "ITScheduler010",
                                    FatalErrorInArgument,
                                    exceptionDescription);
                    }

#ifdef G4VERBOSE
                    if (fVerbose)
                    {
                        G4cout << GetIT(secondary)->GetName() << " ("
                               << secondary->GetTrackID() << ")";
                    }
#endif
                }
            }
            else
            {
#ifdef G4VERBOSE
                if (fVerbose)
                {
                    G4cout << "No product";
                }
#endif
            }
#ifdef G4VERBOSE
            if (fVerbose)
            {
                G4cout << G4endl;
            }
#endif
            if (pTrackA->GetTrackID() == 0 || pTrackB->GetTrackID() == 0)
            {
                G4Track* pTrack = nullptr;
                if (pTrackA->GetTrackID() == 0)
                {
                    pTrack = pTrackA;
                }
                else
                {
                    pTrack = pTrackB;
                }

                G4ExceptionDescription exceptionDescription;
                exceptionDescription
                    << "The problem was found for the reaction between tracks :"
                    << pTrackA->GetParticleDefinition()->GetParticleName() << " ("
                    << pTrackA->GetTrackID() << ") & "
                    << pTrackB->GetParticleDefinition()->GetParticleName() << " ("
                    << pTrackB->GetTrackID() << "). \n";

                if (pTrack->GetStep() == nullptr)
                {
                    exceptionDescription << "Also no step was found"
                                         << " ie track->GetStep() == 0 \n";
                }

                exceptionDescription << "Parent ID of trackA : "
                                     << pTrackA->GetParentID() << "\n";
                exceptionDescription << "Parent ID of trackB : "
                                     << pTrackB->GetParentID() << "\n";

                exceptionDescription
                    << "The ID of one of the reaction track was not setup.";
                G4Exception("G4Scheduler::ComputeInteractionBetweenTracks",
                            "ITScheduler011",
                            FatalErrorInArgument,
                            exceptionDescription);
            }

            if (pChanges->WereParentsKilled())
            {
                pTrackA->SetTrackStatus(fStopAndKill);
                pTrackB->SetTrackStatus(fStopAndKill);

                fpTrackingManager->EndTracking(pTrackA);
                fpTrackingManager->EndTracking(pTrackB);
            }

            pChanges.reset(nullptr);
        }

        fReactionInfo.clear();
    }

    fReactionSet->CleanAllReaction();

    fpTrackContainer->MergeSecondariesWithMainList();
    fpTrackContainer->KillTracks();
}

//______________________________________________________________________________
void G4ITModelProcessor::FindReaction(G4ITReactionSet* pReactionSet,
                                      const double currentStepTime,
                                      const double /*previousStepTime*/,
                                      const bool reachedUserStepTimeLimit)
{
    if (pReactionSet == nullptr || fActiveModels.empty())
    {
        return;
    }

    G4ITReactionPerTrackMap& reactionPerTrackMap = pReactionSet->GetReactionMap();
    G4VITReactionProcess* pReactionProcess = fpActiveModelWithMinTimeStep->GetReactionProcess();

    for (auto tracks_i = reactionPerTrackMap.begin();
         tracks_i != reactionPerTrackMap.end();
         tracks_i = reactionPerTrackMap.begin())
    {
        G4Track* pTrackA = tracks_i->first;
        if (pTrackA->GetTrackStatus() == fStopAndKill)
        {
            continue;
        }

        G4ITReactionPerTrackPtr reactionPerTrack = tracks_i->second;
        G4ITReactionList& reactionList = reactionPerTrack->GetReactionList();

        assert(reactionList.begin() != reactionList.end());

        for (auto it = reactionList.begin(); it != reactionList.end(); it = reactionList.begin())
        {
            G4ITReactionPtr reaction(*it);
            G4Track* pTrackB = reaction->GetReactant(pTrackA);
            if (pTrackB->GetTrackStatus() == fStopAndKill)
            {
                continue;
            }

            if (pTrackB == pTrackA)
            {
                G4ExceptionDescription exceptionDescription;
                exceptionDescription
                    << "The IT reaction process sent back a reaction between trackA and trackB. ";
                exceptionDescription << "The problem is trackA == trackB";
                G4Exception("G4ITModelProcessor::FindReaction",
                            "ITModelProcessor005",
                            FatalErrorInArgument,
                            exceptionDescription);
            }

            pReactionSet->SelectThisReaction(reaction);

            if (pReactionProcess && pReactionProcess->TestReactibility(*pTrackA,
                                                                       *pTrackB,
                                                                       currentStepTime,
                                                                       reachedUserStepTimeLimit))
            {
                auto pReactionChange = pReactionProcess->MakeReaction(*pTrackA, *pTrackB);

                if (pReactionChange)
                {
                    fReactionInfo.push_back(std::move(pReactionChange));
                    break;
                }
            }
        }
    }

    //assert(G4ITReaction::gAll->empty() == true);
}

void G4ITModelProcessor::SetTrack(const G4Track* track)
{
    fpTrack = track;
}

void G4ITModelProcessor::SetModelHandler(G4ITModelHandler* pModelHandler)
{
    if (fInitialized)
    {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription
            << "You are trying to set a new model while the model processor has alreaday be initialized";
        G4Exception("G4ITModelProcessor::SetModelHandler", "ITModelProcessor001",
                    FatalErrorInArgument, exceptionDescription);
    }
    fpModelHandler = pModelHandler;
}

void G4ITModelProcessor::CleanProcessor()
{
    fpTrack = nullptr;
}

bool G4ITModelProcessor::GetComputeTimeStep() const
{
    return fComputeTimeStep;
}

const G4Track* G4ITModelProcessor::GetTrack() const
{
    return fpTrack;
}

void G4ITModelProcessor::SetTrackingManager(G4ITTrackingManager* pTrackingManager)
{
    fpTrackingManager = pTrackingManager;
}


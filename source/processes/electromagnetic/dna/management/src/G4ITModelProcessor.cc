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
// $Id: G4ITModelProcessor.cc 94010 2015-11-05 10:08:33Z gcosmo $
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
//#include "G4ITTimeStepper.hh"
#include "G4ITReaction.hh"
#include "G4ITTrackHolder.hh"
#include "G4ITTrackingManager.hh"
#include "G4UserTimeStepAction.hh"
#include "G4UnitsTable.hh"
#include <vector>

using namespace std;

//#define DEBUG_MEM

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
#endif

G4ITModelProcessor::G4ITModelProcessor()
{
  fpTrack = 0;
  fpModel = 0;
  fInitialized = false;
  fpModelManager = 0;
  fCurrentModel.assign(G4ITType::size(), std::vector<G4VITStepModel*>());

  for(int i = 0; i < (int) G4ITType::size(); i++)
  {
    fCurrentModel[i].assign(G4ITType::size(), 0);
  }
  fUserMinTimeStep = -1.;
  fTSTimeStep = DBL_MAX;
  fpTrackingManager = 0;
  fReactionSet = 0;
  fpTrackContainer = 0;
  fpModelHandler = 0;
  fComputeTimeStep = false;
  fComputeReaction = false;
}

G4ITModelProcessor::~G4ITModelProcessor()
{
  //dtor
  //    if(fpModelHandler) delete fpModelHandler; deleted by G4Scheduler
  fCurrentModel.clear();
  fReactionInfo.clear();
}

// Should not be used
G4ITModelProcessor::G4ITModelProcessor(const G4ITModelProcessor& /*other*/)
{
  //copy ctorr
  fpTrack = 0;
  fpModelHandler = 0;
  fpModel = 0;
  fInitialized = false;
  fpModelManager = 0;
  fUserMinTimeStep = -1.;
  fTSTimeStep = DBL_MAX;
  fpTrackingManager = 0;
  fReactionSet = 0;
  fpTrackContainer = 0;
  fComputeTimeStep = false;
  fComputeReaction = false;
}

void G4ITModelProcessor::RegisterModel(double time, G4VITStepModel* model)
{
  fpModelHandler->RegisterModel(model, time);
}

// Should not be used
G4ITModelProcessor& G4ITModelProcessor::operator=(const G4ITModelProcessor& rhs)
{
  if(this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}
//______________________________________________________________________________
void G4ITModelProcessor::Initialize()
{
  fpModelHandler->Initialize();
  fReactionSet = G4ITReactionSet::Instance();
  fpTrackContainer = G4ITTrackHolder::Instance();
  fInitialized = true;
  fComputeTimeStep = false;
  fComputeReaction = false;
  if(fpModelHandler->GetTimeStepComputerFlag()) fComputeTimeStep = true;
  if(fpModelHandler->GetReactionProcessFlag()) fComputeReaction = true;
}

//_________________________________________________________________________

G4double G4ITModelProcessor::CalculateMinTimeStep(G4double currentGlobalTime,
                                                  G4double definedMinTimeStep)
{

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  MemStat mem_first, mem_second, mem_diff;
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_first = MemoryUsage();
#endif

  fTSTimeStep = DBL_MAX;

  InitializeStepper(currentGlobalTime, definedMinTimeStep);
  // TODO
  // fpMasterModelProcessor -> InitializeStepper(fGlobalTime, fDefinedMinTimeStep) ;

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "\t || MEM || G4Scheduler::CalculateMinTimeStep || After "
  "computing fpModelProcessor -> InitializeStepper, diff is : "
  << mem_diff
  << G4endl;
#endif

//  G4TrackList::iterator fpMainList_i = fpMainList->begin();

  G4TrackManyList* mainList = fpTrackContainer->GetMainList();
  G4TrackManyList::iterator it = mainList->begin();
  G4TrackManyList::iterator end = mainList->end();

  for (; it != end; ++it)
  {
    G4Track * track = *it;

    if (track == 0)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "No track found.";
      G4Exception("G4Scheduler::CalculateMinStep", "ITScheduler006",
          FatalErrorInArgument, exceptionDescription);
      return 0; // makes coverity happy
    }

#ifdef DEBUG
    G4cout << "*_* " << GetIT(track)->GetName()
    << " ID: " << track->GetTrackID()
    << " at time : " << track->GetGlobalTime()
    << G4endl;
#endif

    G4TrackStatus trackStatus = track->GetTrackStatus();
    if (trackStatus == fStopAndKill || trackStatus == fStopButAlive)
    {
      continue;
    }

    // This Extract... was thought for MT mode at the track level
    //ExtractTimeStepperData(fpModelProcessor) ;

    CalculateTimeStep(track, definedMinTimeStep);

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
  if(fpTrack == 0)
  {
    CleanProcessor();
    return;
  }

  const std::vector<std::vector<G4VITStepModel*> >* model = GetCurrentModel();

  for(unsigned i = 0; i < model->size(); ++i)
  {
    for(unsigned j = 0; j < (*model)[i].size(); ++j)
    {
      G4VITStepModel* mod = (*model)[i][j];

      if(mod == 0)
      {
        continue;
      }

      G4VITTimeStepComputer* stepper(mod->GetTimeStepper());

      G4double sampledMinTimeStep(stepper->GetSampledMinTimeStep());
      G4TrackVectorHandle reactants(stepper->GetReactants());

      if(sampledMinTimeStep < fTSTimeStep)
      {
        /*
         // DEBUG SPECIAL CASE
         if(!reactants)
         {
         G4ExceptionDescription exceptionDescription ;
         CalculateMinTimeStep(); // => at least N (N = nb of tracks) loops   exceptionDescription << "No reactants were found by the time stepper.";
         G4Exception("G4Scheduler::ExtractTimeStepperData","ITScheduler007",
         FatalErrorInArgument,exceptionDescription);
         continue ;
         }
         */

        fTSTimeStep = sampledMinTimeStep;
        //fReactingTracks.clear();

        fReactionSet->CleanAllReaction();
        if(bool(reactants))
        {
          // G4cout << "*** (1) G4Scheduler::ExtractTimeStepperData insert
          //          reactants for " << GetIT(track)->GetName() << G4endl;
          // G4cout << "bool(reactants) = " << bool(reactants) << G4endl;
          // G4cout << reactants->size() << G4endl;
          // G4cout << GetIT(reactants->operator[](0))->GetName() << G4endl;

          //  fReactingTracks.insert(make_pair(track, reactants));
          fReactionSet->AddReactions(fTSTimeStep,
                                     const_cast<G4Track*>(fpTrack),
                                     reactants);
          stepper->ResetReactants();
        }
      }
      else if(fTSTimeStep == sampledMinTimeStep)
      {
        /*
         // DEBUG SPECIAL CASE
         if(!reactants)
         {
         G4ExceptionDescription exceptionDescription ;
         exceptionDescription << "No reactants were found by the time stepper.";
         G4Exception("G4Scheduler::ExtractTimeStepperData","ITScheduler008",
         FatalErrorInArgument,exceptionDescription);
         continue ;
         }
         */
        if(bool(reactants))
        {
          // G4cout << "*** (2) G4Scheduler::ExtractTimeStepperData insert
          //          reactants for " << GetIT(track)->GetName() << G4endl;
          // G4cout << "bool(reactants) = " << bool(reactants) << G4endl;
          // G4cout << "trackA : " << GetIT(track)->GetName()
          //        << " ("<< track->GetTrackID() << ")" << G4endl;
          // G4cout << reactants->size() << G4endl;
          // G4cout << GetIT(reactants->operator[](0))->GetName() << G4endl;

          // fReactingTracks.insert(make_pair(track, reactants));

          fReactionSet->AddReactions(fTSTimeStep,
                                     const_cast<G4Track*>(fpTrack),
                                     reactants);
          stepper->ResetReactants();
        }
      }
      else
      {
        if(bool(reactants))
        {
          stepper->ResetReactants();
        }
      }
    }
  }

  CleanProcessor();
}

//______________________________________________________________________________

void G4ITModelProcessor::InitializeStepper(G4double currentGlobalTime,
                                           G4double userMinTime)
{
  // G4cout << "G4ITModelProcessor::InitializeStepper" << G4endl;
  if(fpModelHandler == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "No G4ITModelHandler was passed to the modelProcessor.";
    G4Exception("G4ITModelProcessor::InitializeStepper",
                "ITModelProcessor002",
                FatalErrorInArgument,
                exceptionDescription);
  }
  const std::vector<std::vector<G4ITModelManager*> >* modelManager =
      fpModelHandler->GetAllModelManager();

  if(modelManager == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "No G4ITModelManager was register to G4ITModelHandler.";
    G4Exception("G4ITModelProcessor::InitializeStepper",
                "ITModelProcessor003",
                FatalErrorInArgument,
                exceptionDescription);
  }

  int nbModels1 = modelManager->size();

  G4VITTimeStepComputer::SetTimes(currentGlobalTime, userMinTime);

  // TODO !!!
  //    if( nbModels1 != 1 || (nbModels1 == 1 && !fpModelManager) )
  {
    int nbModels2 = -1;
    G4VITStepModel* model = 0;
    G4ITModelManager* modman = 0;

    for(int i = 0; i < nbModels1; i++)
    {
      nbModels2 = (*modelManager)[i].size();

      for(int j = 0; j <= i; j++)
      {
        modman = (*modelManager)[i][j];

        if(modman == 0) continue;

        model = modman->GetModel(currentGlobalTime);
        G4VITTimeStepComputer* stepper = model->GetTimeStepper();

#if defined (DEBUG_MEM)
        MemStat mem_first, mem_second, mem_diff;
#endif

#if defined (DEBUG_MEM)
        mem_first = MemoryUsage();
#endif

        //      stepper->PrepareForAllProcessors() ;
        stepper->Prepare();

#if defined (DEBUG_MEM)
        mem_second = MemoryUsage();
        mem_diff = mem_second-mem_first;
        G4cout << "\t || MEM || G4ITModelProcessor::InitializeStepper || After computing stepper -> Prepare(), diff is : " << mem_diff << G4endl;
#endif
        fCurrentModel[i][j] = model;
      }
    }

    if(nbModels1 == 1 && nbModels2 == 1)
    {
      fpModelManager = modman;
      fpModel = model;
    }
    else fpModel = 0;
  }
}

//______________________________________________________________________________
void G4ITModelProcessor::CalculateTimeStep(const G4Track* track,
                                           const G4double userMinTimeStep)
{
  // G4cout  << "G4ITModelProcessor::CalculateStep" << G4endl;
  CleanProcessor();
  if(track == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No track was passed to the method (track == 0).";
    G4Exception("G4ITModelProcessor::CalculateStep",
                "ITModelProcessor004",
                FatalErrorInArgument,
                exceptionDescription);
  }
  SetTrack(track);
  fUserMinTimeStep = userMinTimeStep;

  DoCalculateStep();
}

//______________________________________________________________________________
void G4ITModelProcessor::DoCalculateStep()
{
  if(fpModel) // ie only one model has been declared and will be used
  {
    fpModel->GetTimeStepper()->CalculateStep(*fpTrack, fUserMinTimeStep);
  }
  else // ie many models have been declared and will be used
  {
    std::vector<G4VITStepModel*>& model = fCurrentModel[GetIT(fpTrack)
        ->GetITType()];

    for(int i = 0; i < (int) model.size(); i++)
    {
      if(model[i] == 0) continue;
      model[i]->GetTimeStepper()->CalculateStep(*fpTrack, fUserMinTimeStep);
    }
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
  if(fReactionSet->Empty())
  {
    return;
  }

  if(fITStepStatus == eCollisionBetweenTracks)
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

    std::vector<G4ITReactionChange*> * reactionInfo_v = GetReactionInfo();
    std::vector<G4ITReactionChange*>::iterator reactionInfo_i = reactionInfo_v
        ->begin();

    for(; reactionInfo_i != reactionInfo_v->end(); ++reactionInfo_i)
    {
      G4ITReactionChange* changes = (*reactionInfo_i);
      G4Track* trackA = const_cast<G4Track*>(changes->GetTrackA());
      G4Track* trackB = const_cast<G4Track*>(changes->GetTrackB());

      if(trackA == 0 || trackB == 0 || trackA->GetTrackStatus() == fStopAndKill
         || trackB->GetTrackStatus() == fStopAndKill) continue;

      G4int nbSecondaries = changes->GetNumberOfSecondaries();
      const vector<G4Track*>* productsVector = changes->GetfSecondary();

      if(fpUserTimeStepAction)
      {
        fpUserTimeStepAction->UserReactionAction(*trackA,
                                                 *trackB,
                                                 productsVector);
      }

#ifdef G4VERBOSE
      if(fVerbose)
      {
        G4cout << "At time : " << setw(7) << G4BestUnit(fGlobalTime, "Time")
        << " Reaction : " << GetIT(trackA)->GetName() << " ("
        << trackA->GetTrackID() << ") + " << GetIT(trackB)->GetName() << " ("
        << trackB->GetTrackID() << ") -> ";
      }
#endif

      if(nbSecondaries > 0)
      {
        for(int i = 0; i < nbSecondaries; ++i)
        {
#ifdef G4VERBOSE
          if(fVerbose && i != 0) G4cout << " + ";
#endif

          G4Track* secondary = (*productsVector)[i]; //changes->GetSecondary(i);
          fpTrackContainer->_PushTrack(secondary);
          GetIT(secondary)->SetParentID(trackA->GetTrackID(),
                                        trackB->GetTrackID());

          if(secondary->GetGlobalTime() - fGlobalTime > fTimeTolerance)
          {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription << "The time of the secondary should not be bigger than the"
            " current global time."
            << " This may cause synchronization problem. If the process you"
            " are using required "
            << "such feature please contact the developpers." << G4endl<< "The global time in the step manager : "
            << G4BestUnit(fGlobalTime,"Time")
            << G4endl
            << "The global time of the track : "
            << G4BestUnit(secondary->GetGlobalTime(),"Time")
            << G4endl;

            G4Exception("G4Scheduler::ComputeInteractionBetweenTracks",
                        "ITScheduler010",
                        FatalErrorInArgument,
                        exceptionDescription);
          }

#ifdef G4VERBOSE
          if(fVerbose)
            G4cout << GetIT(secondary)->GetName() << " ("
                   << secondary->GetTrackID() << ")";
#endif
        }
      }
      else
      {
#ifdef G4VERBOSE
        if(fVerbose) G4cout << "No product";
#endif
      }
#ifdef G4VERBOSE
      if(fVerbose) G4cout << G4endl;
#endif
      if(trackA->GetTrackID() == 0 || trackB->GetTrackID() == 0)
      {
        G4Track* track = 0;
        if(trackA->GetTrackID() == 0) track = trackA;
        else track = trackB;

        G4ExceptionDescription exceptionDescription;
        exceptionDescription
            << "The problem was found for the reaction between tracks :"
            << trackA->GetParticleDefinition()->GetParticleName() << " ("
            << trackA->GetTrackID() << ") & "
            << trackB->GetParticleDefinition()->GetParticleName() << " ("
            << trackB->GetTrackID() << "). \n";

        if(track->GetStep() == 0)
        {
          exceptionDescription << "Also no step was found"
                               << " ie track->GetStep() == 0 \n";
        }

        exceptionDescription << "Parent ID of trackA : "
                             << trackA->GetParentID() << "\n";
        exceptionDescription << "Parent ID of trackB : "
                             << trackB->GetParentID() << "\n";

        exceptionDescription
            << "The ID of one of the reaction track was not setup.";
        G4Exception("G4Scheduler::ComputeInteractionBetweenTracks",
                    "ITScheduler011",
                    FatalErrorInArgument,
                    exceptionDescription);
      }

      if(changes->WereParentsKilled())
      {
        // DEBUG
        // G4cout << "Erasing tracks : "
        //  << "trackA at time : " << G4BestUnit(trackA->GetGlobalTime() , "Time")
        //  << "\t trackB at time : "<< G4BestUnit(trackB->GetGlobalTime(), "Time")
        //  << "\t GlobalTime : " << G4BestUnit(fGlobalTime, "Time")
        //  << G4endl;

        trackA->SetTrackStatus(fStopAndKill);
        trackB->SetTrackStatus(fStopAndKill);

//        G4TrackList::Pop(trackA);
//        G4TrackList::Pop(trackB);
        fpTrackingManager->EndTracking(trackA);
        fpTrackingManager->EndTracking(trackB);
      }

      delete changes;
    }

    reactionInfo_v->clear(); // never null because point to a non-pointer member of ModelProcessor
  }
  // DEBUG
  //else
  //{
  //    G4cout << "fInteractionStep == true" << G4endl ;
  //}

//  fReactingTracks.clear();
  fReactionSet->CleanAllReaction();

  fpTrackContainer->MergeSecondariesWithMainList();
  fpTrackContainer->KillTracks();
}

//______________________________________________________________________________
void G4ITModelProcessor::FindReaction(G4ITReactionSet* reactionSet,
                                      const double currentStepTime,
                                      const double previousStepTime,
                                      const bool reachedUserStepTimeLimit)
{
  // DEBUG
  //    G4cout << "G4ITReactionManager::FindReaction" << G4endl;
  //if (tracks == 0) return;
  if(reactionSet == 0) return;
  if(fpModelHandler->GetAllModelManager()->empty()) return;

  G4ITReactionPerTrackMap& reactionPerTrackMap = reactionSet->GetReactionMap();

  std::map<G4Track*, G4ITReactionPerTrackPtr, compTrackPerID>::iterator tracks_i =
      reactionPerTrackMap.begin();

  //std::map<G4Track*, G4TrackVectorHandle, compTrackPerID>::iterator tracks_i = tracks->begin();

  //    G4cout << "G4ITModelProcessor::FindReaction at step :" << G4ITTimeStepper::Instance()->GetNbSteps() << G4endl;

  //  for (tracks_i = tracks->begin(); tracks_i != tracks->end(); tracks_i++)
  for(tracks_i = reactionPerTrackMap.begin();
      tracks_i != reactionPerTrackMap.end();
      tracks_i = reactionPerTrackMap.begin())
  {
    //G4cout << "here" << G4endl;
    G4Track* trackA = tracks_i->first;
    if(trackA->GetTrackStatus() == fStopAndKill)
    {
      //G4cout << "continue 1" << G4endl;
      continue;
    }

    G4ITReactionPerTrackPtr reactionPerTrack = tracks_i->second;
    G4ITReactionList& reactionList = reactionPerTrack->GetReactionList();

    G4IT* ITA = GetIT(trackA);
    G4ITType ITypeA = ITA->GetITType();

    const std::vector<G4VITStepModel*> model = fCurrentModel[ITypeA];

    G4Track* trackB = 0;
    G4ITType ITypeB(-1);
    G4VITReactionProcess* process = 0;
    G4ITReactionChange* changes = 0;

    assert(reactionList.begin() != reactionList.end());

    for(G4ITReactionList::iterator it = reactionList.begin();
        it != reactionList.end(); it = reactionList.begin())
    {
      G4ITReactionPtr reaction(*it);
      trackB = reaction->GetReactant(trackA);
      if(trackB->GetTrackStatus() == fStopAndKill)
      {
        //G4cout << "continue 2" << G4endl;
        continue;
      }

      if(trackB == trackA)
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

      G4IT* ITB = GetIT(trackB);
      G4ITType ITypeBtmp = ITB->GetITType();

      if(ITypeB != ITypeBtmp)
      {
        ITypeB = ITypeBtmp;

        if(model[ITypeB]) process = model[ITypeB]->GetReactionProcess();
      }

      reactionSet->SelectThisReaction(reaction);

      if(process && process->TestReactibility(*trackA,
                                              *trackB,
                                              currentStepTime,
                                              previousStepTime,
                                              reachedUserStepTimeLimit))
      {
        changes = process->MakeReaction(*trackA, *trackB);
      }

      if(changes)
      {
        fReactionInfo.push_back(changes);

        //        G4cout << "pushing reaction for trackA (" << trackA->GetTrackID() << ") and trackB ("
        //             << trackB->GetTrackID() << ")" << G4endl;
        //
        //        G4cout << "nb of secondaries : " << changes->GetNumberOfSecondaries() << G4endl;
        //        G4cout << "with track 0 = " << changes->GetSecondary(0) << G4endl;

        process->ResetChanges();
        changes = 0;

        break;
      }
    }
  }

  //assert(G4ITReaction::gAll->empty() == true);
}

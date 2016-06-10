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
// $Id: G4ITModelProcessor.cc 90769 2015-06-09 10:33:41Z gcosmo $
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

//#define DEBUG_MEM

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
#endif

G4ThreadLocal std::map<const G4Track*, G4bool> *G4ITModelProcessor::fHasReacted =
    0;

G4ITModelProcessor::G4ITModelProcessor()
{
  //ctor
  if (!fHasReacted) fHasReacted = new std::map<const G4Track*, G4bool>;
  fpTrack = 0;
  fpModelHandler = 0;
  fpModel = 0;
  fInitialized = false;
  fpModelManager = 0;
  fCurrentModel.assign(G4ITType::size(), std::vector<G4VITStepModel*>());

  for (int i = 0; i < (int) G4ITType::size(); i++)
  {
    fCurrentModel[i].assign(G4ITType::size(), 0);
  }
  fUserMinTimeStep = -1.;
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
}

// Should not be used
G4ITModelProcessor& G4ITModelProcessor::operator=(const G4ITModelProcessor& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}
//______________________________________________________________________________
void G4ITModelProcessor::Initialize()
{
  fpModelHandler->Initialize();
  fInitialized = true;
}

//______________________________________________________________________________
void G4ITModelProcessor::InitializeStepper(const G4double& currentGlobalTime,
                                           const G4double& userMinTime)
{
  // G4cout << "G4ITModelProcessor::InitializeStepper" << G4endl;
  if (fpModelHandler == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
    << "No G4ITModelHandler was passed to the modelProcessor.";
    G4Exception("G4ITModelProcessor::InitializeStepper", "ITModelProcessor002",
                FatalErrorInArgument, exceptionDescription);
  }
  const std::vector<std::vector<G4ITModelManager*> >* modelManager =
      fpModelHandler->GetAllModelManager();

  if (modelManager == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
    << "No G4ITModelManager was register to G4ITModelHandler.";
    G4Exception("G4ITModelProcessor::InitializeStepper", "ITModelProcessor003",
                FatalErrorInArgument, exceptionDescription);
  }

  int nbModels1 = modelManager->size();

  G4VITTimeStepComputer::SetTimes(currentGlobalTime, userMinTime);

  // TODO !!!
  //    if( nbModels1 != 1 || (nbModels1 == 1 && !fpModelManager) )
  {
    int nbModels2 = -1;
    G4VITStepModel* model = 0;
    G4ITModelManager* modman = 0;

    for (int i = 0; i < nbModels1; i++)
    {
      nbModels2 = (*modelManager)[i].size();

      for (int j = 0; j <= i; j++)
      {
        modman = (*modelManager)[i][j];

        if (modman == 0) continue;

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

    if (nbModels1 == 1 && nbModels2 == 1)
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
  if (track == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No track was passed to the method (track == 0).";
    G4Exception("G4ITModelProcessor::CalculateStep", "ITModelProcessor004",
                FatalErrorInArgument, exceptionDescription);
  }
  SetTrack(track);
  fUserMinTimeStep = userMinTimeStep;

  DoCalculateStep();
}

//______________________________________________________________________________
void G4ITModelProcessor::DoCalculateStep()
{
  if (fpModel) // ie only one model has been declared and will be used
  {
    fpModel->GetTimeStepper()->CalculateStep(*fpTrack, fUserMinTimeStep);
  }
  else // ie many models have been declared and will be used
  {
    std::vector<G4VITStepModel*>& model = fCurrentModel[GetIT(fpTrack)
                                                        ->GetITType()];

    for (int i = 0; i < (int) model.size(); i++)
    {
      if (model[i] == 0) continue;
      model[i]->GetTimeStepper()->CalculateStep(*fpTrack, fUserMinTimeStep);
    }
  }
}

//______________________________________________________________________________
void G4ITModelProcessor::FindReaction(
    G4ITReactionSet* reactionSet,
    const double currentStepTime,
    const double previousStepTime,
    const bool reachedUserStepTimeLimit)
{
  // DEBUG
  //    G4cout << "G4ITReactionManager::FindReaction" << G4endl;
  //if (tracks == 0) return;
  if(reactionSet == 0) return;
  if (fpModelHandler->GetAllModelManager()->empty()) return;

  G4ITReactionPerTrackMap& reactionPerTrackMap = reactionSet->GetReactionMap();

  std::map<G4Track*, G4ITReactionPerTrackPtr, compTrackPerID>::iterator tracks_i = reactionPerTrackMap.begin();

  //std::map<G4Track*, G4TrackVectorHandle, compTrackPerID>::iterator tracks_i = tracks->begin();

  //    G4cout << "G4ITModelProcessor::FindReaction at step :" << G4ITTimeStepper::Instance()->GetNbSteps() << G4endl;

  //  for (tracks_i = tracks->begin(); tracks_i != tracks->end(); tracks_i++)
  for (tracks_i = reactionPerTrackMap.begin();
       tracks_i != reactionPerTrackMap.end() ;
       tracks_i = reactionPerTrackMap.begin())
  {
    //G4cout << "here" << G4endl;
    G4Track* trackA = tracks_i->first;
    if (trackA->GetTrackStatus() == fStopAndKill)
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

    for(G4ITReactionList::iterator it = reactionList.begin() ;
        it != reactionList.end() ; it = reactionList.begin() )
    {
      G4ITReactionPtr reaction(*it);
      trackB = reaction->GetReactant(trackA);
      if(trackB->GetTrackStatus() == fStopAndKill)
      {
        //G4cout << "continue 2" << G4endl;
        continue;
      }

      if (trackB == trackA)
      {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription
        << "The IT reaction process sent back a reaction between trackA and trackB. ";
        exceptionDescription << "The problem is trackA == trackB";
        G4Exception("G4ITModelProcessor::FindReaction", "ITModelProcessor005",
                    FatalErrorInArgument, exceptionDescription);
      }

      G4IT* ITB = GetIT(trackB);
      G4ITType ITypeBtmp = ITB->GetITType();

      if (ITypeB != ITypeBtmp)
      {
        ITypeB = ITypeBtmp;

        if (model[ITypeB]) process = model[ITypeB]->GetReactionProcess();
      }
      
      reactionSet->SelectThisReaction(reaction);

      if (process && process->TestReactibility(*trackA, *trackB,
                                               currentStepTime,
                                               previousStepTime,
                                               reachedUserStepTimeLimit))
      {
        changes = process->MakeReaction(*trackA, *trackB);
      }

      if (changes)
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

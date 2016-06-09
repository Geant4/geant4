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

#include "G4ITStepManager.hh"
#include "G4ITModelProcessor.hh"
#include "G4ITStepProcessor.hh"
#include "G4IT.hh"
#include "G4ITReactionChange.hh"
#include "G4ITModelHandler.hh"
#include "G4VITModel.hh"
#include "G4UserReactionAction.hh"
#include "G4ITTrackingManager.hh"
#include "G4TrackingInformation.hh"
#include "G4UnitsTable.hh"

using namespace std;

auto_ptr<G4ITStepManager> G4ITStepManager::fgStepManager(0) ;

template<typename T>
inline bool IsInf(T value)
{
    return std::numeric_limits<T>::has_infinity &&
            value == std::numeric_limits<T>::infinity();
}
//_________________________________________________________________________

G4ITStepManager* G4ITStepManager::Instance()
{
    if(fgStepManager.get() == 0)
        fgStepManager = auto_ptr<G4ITStepManager>(new G4ITStepManager());
    return fgStepManager.get() ;
} 
//_________________________________________________________________________

void G4ITStepManager::DeleteInstance()
{
    if(fgStepManager.get()) fgStepManager.reset();
}
//_________________________________________________________________________

G4ITStepManager::G4ITStepManager()
{
    fpMainList = 0 ;
    fpWaitingList = 0 ;

    fpUserTimeSteps = 0;

    fTimeStep = DBL_MAX ;
    fTSTimeStep = DBL_MAX;
    fILTimeStep = DBL_MAX;
    fPreviousStepTime = DBL_MAX;

    fComputeTimeStep = false;
    fComputeReaction = false;
    fTimeTolerance = 1 * picosecond;
    fEndTime = 1 * microsecond ;
    fGlobalTime = -1 ;
    fInteractionStep = true ;
    fUsePreDefinedTimeSteps = false;

    fpStepProcessor = 0 ;
    fpMasterStepProcessor = 0 ;

    fpModelProcessor= 0;
    fpMasterModelProcessor= 0;

    fNbSteps = 0;

    fRunning = false;
    fInitialized = false;

    fpUserReactionAction = 0;
    fpTrackingManager = 0;

    fNbTracks = -1;
    fpModelHandler = new G4ITModelHandler();
    fpTrackingManager = new G4ITTrackingManager();

    fVerbose = 0;
    fDefinedMinTimeStep = -1.;
    fReachedUserTimeLimit = false;
    fTmpEndTime = -1.;
    fTmpGlobalTime = -1.;
}
//_________________________________________________________________________

G4ITStepManager::~G4ITStepManager()
{
    if(fpMasterStepProcessor) delete fpMasterStepProcessor ;
    if(fpMasterModelProcessor) delete fpMasterModelProcessor ;
    G4ITTypeManager::DeleteInstance();
    ClearList();
    if(fpTrackingManager) delete fpTrackingManager;
    fgStepManager.release();
}
//_________________________________________________________________________

void G4ITStepManager::ClearList()
{
    if(fNbTracks == 0) return;
    if(fpMainList)
    {
        delete fpMainList;
        fpMainList = 0;
    }

    if(fpWaitingList)
    {
        delete fpWaitingList;
        fpWaitingList = 0;
    }

    if(!fDelayedList.empty())
    {
        std::map<double,G4TrackList* >::iterator fDelayedList_i = fDelayedList.begin() ;

        for(; fDelayedList_i != fDelayedList.end() ; fDelayedList_i++)
        {
            if(fDelayedList_i->second)
                delete (fDelayedList_i->second);
            fDelayedList_i->second = 0 ;
        }
        fDelayedList.clear();
    }
    fNbTracks = -1;
}

//_________________________________________________________________________

void G4ITStepManager::Initialize()
{
    if(fpStepProcessor)
    {
        delete fpStepProcessor ;
    }
    if(fpMasterModelProcessor)
    {
        delete fpMasterModelProcessor;
    }

    //______________________________________________________________

    fpMasterModelProcessor = new G4ITModelProcessor();
    fpMasterModelProcessor->SetModelHandler(fpModelHandler);
    fpModelProcessor      = fpMasterModelProcessor;

    //______________________________________________________________

    fpMasterStepProcessor = new G4ITStepProcessor();
    fpMasterStepProcessor->SetTrackingManager(fpTrackingManager);
    fpStepProcessor = fpMasterStepProcessor ;

    if(fpModelHandler->GetTimeStepComputerFlag()) fComputeTimeStep = true;
    if(fpModelHandler->GetReactionProcessFlag())  fComputeReaction = true;
    //______________________________________________________________

    fInitialized = true ;
}
//_________________________________________________________________________

void G4ITStepManager::Reset()
{
    fTimeStep = DBL_MAX ;
    fTSTimeStep = DBL_MAX;
    fILTimeStep = DBL_MAX;
    fPreviousStepTime = DBL_MAX;
    fGlobalTime = -1 ;
    fInteractionStep = true ;

    fNbSteps = 0;
}
//_________________________________________________________________________

void G4ITStepManager::Process()
{

#ifdef G4VERBOSE
    if(fVerbose)
    {
        G4cout << "*** G4ITStepManager starts processing "<< G4endl;
        if(fVerbose > 1)
            G4cout  << "______________________________________________________________________" << G4endl;
    }
#endif

    fpTrackingManager->Initialize();
    fpMasterModelProcessor->Initialize();
    fpMasterStepProcessor->Initialize();

    // ___________________
    fRunning = true ;
    fPreviousStepTime = DBL_MAX;
    if(fpUserReactionAction) fpUserReactionAction->StartProcessing();

    if(fDelayedList.empty() == false)     SynchronizeTracks() ;
    DoProcess() ;

#ifdef G4VERBOSE
    if(fVerbose)
    {
        G4cout << "*** G4ITStepManager ends at time : " << G4BestUnit(fGlobalTime,"Time") << G4endl;
        G4cout  << "___________________________________" << G4endl;
    }
#endif

    if(fpUserReactionAction) fpUserReactionAction->EndProcessing();
    EndTracking();

    // ___________________
    fRunning = false;
}
//_________________________________________________________________________

void G4ITStepManager::SynchronizeTracks()
{
    if(fpWaitingList != 0)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "There is a waiting track list (fpWaitingList != 0).";
        exceptionDescription << " When G4ITStepManager::SynchronizeTracks() is called, ";
        exceptionDescription << "no more tracks should remain in the fpWaitingList.";
        G4Exception("G4ITStepManager::SynchronizeTracks","ITStepManager002",
                    FatalErrorInArgument,exceptionDescription);
    }

    // If ther's not yet the main list, retrieve the max timings
    // of the delayed tracks and consider the list of the correspondant tracks
    // as the main list.

    if(!fpMainList)
    {
        std::map<double,G4TrackList* >::iterator
                fDelayedList_i = fDelayedList.end();
        fDelayedList_i -- ;

        if(fDelayedList_i->second && ! fDelayedList_i->second->empty())
        {
            fGlobalTime = fDelayedList_i->first ;
            fpMainList = fDelayedList_i->second ;

            //            DEBUG
            //            G4cout << "Selected Global time : " << fGlobalTime  << G4endl;
            //            G4cout << fpMainList->front()->GetGlobalTime() << G4endl;

            fDelayedList_i->second = 0 ;
        }
        fDelayedList[fDelayedList_i->first] = 0 ;
        fDelayedList.erase(fDelayedList_i);

        // If There are no more delayed tracks
        if(fDelayedList.empty())
        {
            return ;
        }
    }

    // Backup main list and time feature
    // Reminder : the main list here, should
    // have the biggest global time
    
    fpWaitingList   = fpMainList ;
    fTmpGlobalTime  = fGlobalTime ;
    fTmpEndTime     = fEndTime ;

    std::map<double,G4TrackList* >::iterator fDelayedList_i = fDelayedList.begin();
    std::map<double,G4TrackList* >::iterator nextfDelayedList_i;
    if(fDelayedList_i != fDelayedList.end())
    {
        // Retrieve the first delayed list
        fpMainList      = fDelayedList_i->second ;
        fGlobalTime     = fDelayedList_i->first;

        do
        {
            if((fDelayedList_i->second) == 0)   continue ;
            if(fDelayedList_i->second->empty()) continue ;

            if(fDelayedList_i->first > fTmpGlobalTime)
            {

                G4ExceptionDescription exceptionDescription ;
                exceptionDescription << "The next time for delayed tracks is bigger than the maximum allowed simulation time.";
                exceptionDescription << "Next time to consider = " << G4BestUnit(fDelayedList_i->first,"Time") << "\n";
                exceptionDescription << "G4ITStepManager::fTmpGlobalTime (should be the biggest time in the simulation) = "
                                     << G4BestUnit(fTmpGlobalTime,"Time")<< G4endl;

                G4Exception("G4ITStepManager::SynchronizeTracks","ITStepManager003",
                            FatalErrorInArgument,exceptionDescription);
            }

            nextfDelayedList_i = fDelayedList_i ;
            nextfDelayedList_i++;

            if(nextfDelayedList_i != fDelayedList.end())
            {
                fEndTime = (nextfDelayedList_i)->first ;
            }
            else
            {
                fEndTime = fTmpEndTime ;
            }

            DoProcess();

            if(nextfDelayedList_i != fDelayedList.end())
            {
                // Merge next list with current one
                G4TrackList* nextList = nextfDelayedList_i->second;
                nextList -> transferTo (fpMainList);
                delete nextList ;
                nextfDelayedList_i->second = 0 ;
            }
            else
            {

                // Merge first list put in waitingTracks
                fpWaitingList -> transferTo (fpMainList);
                delete fpWaitingList ;
                fpWaitingList = 0 ;
            }
            fDelayedList[fDelayedList_i->first] = 0 ;

            fDelayedList_i++;
        }
        while (fDelayedList_i != fDelayedList.end());
        fDelayedList.clear();
    }
    else
    {
        // If no delayed list, then reput the main list to its place
        fpMainList = fpWaitingList ;
        fpWaitingList = 0;
    }
    fGlobalTime    = fTmpGlobalTime ;
    fEndTime       = fTmpEndTime ;

}
//_________________________________________________________________________

void G4ITStepManager::DoProcess()
// We split it from the Process() method to avoid repeating code in SynchronizeTracks
{
    while(fGlobalTime <= fEndTime &&  ! (fpMainList -> empty()) )
    {
        Stepping() ;
    }

#ifdef G4VERBOSE
    if(fVerbose > 1)
        G4cout << "*** G4ITStepManager has finished processing a track list at time : " << G4BestUnit(fGlobalTime,"Time") << G4endl;
#endif
}
//_________________________________________________________________________

void G4ITStepManager::Stepping()
{
    fTimeStep = DBL_MAX ;

    fTSTimeStep = DBL_MAX;
    fILTimeStep = DBL_MAX;

    fInteractionStep = false ;
    fReachedUserTimeLimit = false;

    // Start of step
#ifdef G4VERBOSE
    if(fVerbose > 1)
    {
//        G4cout  << "\33[1;31m";
        G4cout  << "*** Start Of Step N°" << fNbSteps+1 << " ***" << G4endl;
        G4cout  << "Current Global time : " << G4BestUnit(fGlobalTime, "Time") <<G4endl;
//        G4cout  << "\33[0m";
    }
#endif

    if(fUsePreDefinedTimeSteps)
    {
        FindUserPreDefinedTimeStep();

#ifdef G4VERBOSE
        if(fVerbose > 1)
        {
//            G4cout  << "\33[1;31m";
            G4cout  << "*** At time : " << G4BestUnit(fGlobalTime, "Time")
                    << " the chosen user time step is : " << G4BestUnit(fDefinedMinTimeStep, "Time")
                    << " ***" << G4endl;
//            G4cout  << "\33[0m";
        }
#endif
    }
    else fDefinedMinTimeStep = 0 ;


    if(fComputeTimeStep)
    {
        CalculateMinTimeStep() ; // => at least N (N = nb of tracks) loops
    }
    else
    {
        fTSTimeStep = fDefinedMinTimeStep ;
    }

#ifdef G4VERBOSE
    if(fVerbose > 1)
    {
//        G4cout  << "\33[1;31m";
        G4cout  << "*** Time stepper returned : " << G4BestUnit(fTSTimeStep, "Time") << " ***" << G4endl;
//        G4cout  << "\33[0m";
    }
#endif


    // Call IL even if fTSTimeStep == 0
    // if fILTimeStep == 0 give the priority to DoIt processes
    ComputeInteractionLength(); // => at least N loops
    // All process returns the physical step of interaction
    // using the approximation : E = cste along the step
    // Only the transportation calculates the corresponding
    // time step

#ifdef G4VERBOSE
    if(fVerbose > 1)
    {
//        G4cout  << "\33[1;31m";
        G4cout  << "*** The minimum time returned by the processes is : " << G4BestUnit(fILTimeStep, "Time") << " ***" << G4endl;
//        G4cout  << "\33[0m";
    }
#endif


    if(fILTimeStep <= fTSTimeStep) // Give the priority to the IL
    {
        fInteractionStep = true ;
        fReactingTracks.clear(); // Give the priority to the IL

        fTimeStep = fILTimeStep;
    }
    else
    {
        fInteractionStep = false ;
        if(fLeadingTracks.empty() == false)
        {
            std::vector<G4Track*>::iterator fLeadingTracks_i = fLeadingTracks.begin();

            while(fLeadingTracks_i != fLeadingTracks.end())
            {
                G4Track* track = *fLeadingTracks_i;
                if(track)
                {
                    G4IT* ITrack = GetIT(*fLeadingTracks_i) ;
                    if(ITrack)
                    {
                        GetIT(*fLeadingTracks_i)->GetTrackingInfo()->SetLeadingStep(false);
                    }
                }

                fLeadingTracks_i++;
                continue ;
            }

            fLeadingTracks.clear();
        }

        fTimeStep = fTSTimeStep;
    }

    fReachedUserTimeLimit = ((fTimeStep <= fDefinedMinTimeStep) ||
                             ((fTimeStep > fDefinedMinTimeStep) &&
                              fabs(fTimeStep-fDefinedMinTimeStep) < fTimeTolerance)) ?
                true : false;

    if(fpUserReactionAction) fpUserReactionAction->TimeStepAction();

    // if fTSTimeStep > 0 => needs to call the transportation process
    // if fILTimeStep < fTSTimeStep => call DoIt processes
    // if fILTimeStep == fTSTimeStep => give the priority to the DoIt processes
    if(fTSTimeStep > 0 || fILTimeStep <= fTSTimeStep)
    {
        DoIt();

        MergeSecondariesWithMainList();
        KillTracks();
    }

    fGlobalTime += fTimeStep ;

    ComputeTrackReaction();

    MergeSecondariesWithMainList();

    KillTracks();

    fNbSteps++;

    fPreviousStepTime = fTimeStep;

    // End of step
#ifdef G4VERBOSE
    if(fVerbose > 1)
    {
//        G4cout  << "\33[1;31m";
        G4cout  << "*** End Of Step N°" << fNbSteps << " ***" << G4endl;
        G4cout  << "Previous global time : " << G4BestUnit(fGlobalTime- fTimeStep, "Time")
                << "\t Chosen step time :"
                << G4BestUnit(fTimeStep,"Time")<<G4endl;
        if(fReachedUserTimeLimit) G4cout << "It has also reached the user time limit" << G4endl;
        G4cout  << "Next global time : " << G4BestUnit(fGlobalTime,"Time") << G4endl;
        G4cout  << "______________________________________________________________________"<< G4endl;
//        G4cout  << "\33[0m" ;
    }
#endif

}
//_________________________________________________________________________

void G4ITStepManager::FindUserPreDefinedTimeStep()
{

    if(fpUserTimeSteps == 0)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "You are asking to use user defined steps but you did not give any.";
        G4Exception("G4ITStepManager::FindUserPreDefinedTimeStep","ITStepManager004",
                    FatalErrorInArgument,exceptionDescription);
        return ; // makes coverity happy
    }
    map<double, double>::iterator fpUserTimeSteps_i   = fpUserTimeSteps->upper_bound(fGlobalTime) ;
    map<double, double>::iterator fpUserTimeSteps_low = fpUserTimeSteps->lower_bound(fGlobalTime) ;

    // DEBUG
    //        G4cout << "fGlobalTime : " << G4BestUnit(fGlobalTime,"Time") << G4endl;
    //        G4cout << "fpUserTimeSteps_i : "
    //        <<"<"<<G4BestUnit(fpUserTimeSteps_i->first,"Time")<<", "<< G4BestUnit(fpUserTimeSteps_i->second,"Time")<<">"
    //        << "\t fpUserTimeSteps_low : "
    //        <<"<"<<G4BestUnit(fpUserTimeSteps_low->first,"Time")<<", "<< G4BestUnit(fpUserTimeSteps_low->second,"Time")<<">"
    //        << G4endl;


    if(fpUserTimeSteps_i == fpUserTimeSteps->end())
    {
        fpUserTimeSteps_i --;
    }
    else if(fabs(fGlobalTime - fpUserTimeSteps_low->first) < fTimeTolerance )
    {
        // Case : fGlobalTime = X picosecond
        // and fpUserTimeSteps_low->first = X picosecond
        // but the precision is not good enough
        fpUserTimeSteps_i = fpUserTimeSteps_low;
    }
    else if(fpUserTimeSteps_i == fpUserTimeSteps_low)
    {
        // "Normal" cases
        fpUserTimeSteps_i--;
    }
    else
    {
        fpUserTimeSteps_i = fpUserTimeSteps_low;
    }

    fDefinedMinTimeStep = fpUserTimeSteps_i->second ;
}
//_________________________________________________________________________

void G4ITStepManager::CalculateMinTimeStep()
{

    if(fpMasterModelProcessor == 0)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "There is no G4ITModelProcessor to hande IT reaction. ";
        exceptionDescription << "You probably did not initialize the G4ITStepManager. ";
        exceptionDescription << "Just do G4ITStepManager::Instance()->Initialize(); ";
        exceptionDescription << " but only after initializing the run manager.";
        G4Exception("G4ITStepManager::CalculateMinStep","ITStepManager005",
                    FatalErrorInArgument,exceptionDescription);
        return ; // makes coverity happy
    }

    fpMasterModelProcessor -> InitializeStepper(fGlobalTime, fDefinedMinTimeStep) ;

    G4TrackList::iterator      fpMainList_i = fpMainList->begin();
    for( ; fpMainList_i != fpMainList->end() ; fpMainList_i++)
    {
        G4Track * track  = *fpMainList_i ;

        if(track == 0)
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "No track found.";
            G4Exception("G4ITStepManager::CalculateMinStep","ITStepManager006",
                        FatalErrorInArgument,exceptionDescription);
            return ; // makes coverity happy
        }

        G4TrackStatus trackStatus = track->GetTrackStatus();
        if(trackStatus == fStopAndKill || trackStatus == fStopButAlive)
        {
            continue ;
        }
        ExtractTimeStepperData(fpModelProcessor) ;

        fpModelProcessor->CalculateTimeStep(track, fDefinedMinTimeStep);
    }

    ExtractTimeStepperData(fpModelProcessor);
}
//_________________________________________________________________________

void G4ITStepManager::ExtractTimeStepperData(G4ITModelProcessor* MP)
{
    G4Track* track = (G4Track*) MP -> GetTrack() ;
    if(track == 0)
    {
        MP->CleanProcessor();
        return ;
    }

    const std::vector<std::vector<G4VITModel*> >* model = MP->GetCurrentModel();

    for(unsigned i = 0 ; i < model->size() ; i++)
    {
        for(unsigned j = 0 ; j < (*model)[i].size() ; j++)
        {
            G4VITModel* mod = (*model)[i][j];

            if(mod == 0)
            {
                continue;
            }

            G4VITTimeStepper* stepper (mod -> GetTimeStepper()) ;

            G4double sampledMinTimeStep (stepper -> GetSampledMinTimeStep()) ;
            G4TrackVectorHandle reactants (stepper -> GetReactants());

            if(sampledMinTimeStep < fTSTimeStep)
            {
                /*
                // DEBUG SPECIAL CASE
                if(!reactants)
                {
                    G4ExceptionDescription exceptionDescription ;
                    exceptionDescription << "No reactants were found by the time stepper.";
                    G4Exception("G4ITStepManager::ExtractTimeStepperData","ITStepManager007",
                                FatalErrorInArgument,exceptionDescription);
                    continue ;
                }
*/

                fTSTimeStep = sampledMinTimeStep ;
                fReactingTracks.clear();
                if(bool(reactants) )
                {
                    fReactingTracks.insert(make_pair(track, reactants));
                    stepper -> ResetReactants();
                }
            }
            else if(fTSTimeStep == sampledMinTimeStep )
            {
                /*
                // DEBUG SPECIAL CASE
                if(!reactants)
                {
                    G4ExceptionDescription exceptionDescription ;
                    exceptionDescription << "No reactants were found by the time stepper.";
                    G4Exception("G4ITStepManager::ExtractTimeStepperData","ITStepManager008",
                                FatalErrorInArgument,exceptionDescription);
                    continue ;
                }
*/
                if(bool(reactants) )
                {
                    fReactingTracks.insert(make_pair(track, reactants));
                    stepper -> ResetReactants();
                }
            }
            else
            {
                if(bool(reactants) )
                {
                    stepper -> ResetReactants();
                }
            }
        }
    }

    MP->CleanProcessor();
}
//_________________________________________________________________________

void G4ITStepManager::ComputeInteractionLength()
{
    G4TrackList::iterator fpMainList_i = fpMainList->begin();

    for( ; fpMainList_i != fpMainList->end() ; fpMainList_i++)
    {
        G4Track * track  = *fpMainList_i ;

        fpStepProcessor -> DefinePhysicalStepLength(track) ;

        if(fpStepProcessor->GetTrack()) ExtractILData(fpStepProcessor) ;
    }
}
//_________________________________________________________________________

void G4ITStepManager::ExtractILData(G4ITStepProcessor* SP)
{
    G4Track* track = SP->GetTrack();
    if(! track)
    {
        SP->CleanProcessor();
        return ;
    }
    if(track->GetTrackStatus() == fStopAndKill)
    {
        fpMainList->pop(track);
        EndTracking(track);
        return ;
    }

    if(IsInf(SP -> GetInteractionTime()))
    {
        SP->CleanProcessor();
        return;
    }
    else if( SP -> GetInteractionTime() < fILTimeStep)
    {
        if(fLeadingTracks.empty() == false)
        {
            std::vector<G4Track*>::iterator fLeadingTracks_i = fLeadingTracks.begin();

            while(fLeadingTracks_i != fLeadingTracks.end())
            {
                G4Track* track = *fLeadingTracks_i;
                if(track)
                {
                    G4IT* ITrack = GetIT(*fLeadingTracks_i) ;
                    if(ITrack)
                    {
                        GetIT(*fLeadingTracks_i)->GetTrackingInfo()->SetLeadingStep(false);
                    }
                }

                fLeadingTracks_i++;
                continue ;
            }

            fLeadingTracks.clear();
        }

        //        G4cout << "Will set leading step to true for time  :"
        //               << SP -> GetInteractionTime() << " against fTimeStep : "
        //               << fTimeStep << " the trackID is : " << track->GetTrackID() << G4endl;

        fILTimeStep =  SP -> GetInteractionTime() ;

        GetIT(track)->GetTrackingInfo()->SetLeadingStep(true);
        fLeadingTracks.push_back(track);
    }
    else if(fILTimeStep == SP -> GetInteractionTime() )
    {

        //        G4cout << "Will set leading step to true for time  :"
        //               << SP -> GetInteractionTime() << " against fTimeStep : "
        //               << fTimeStep << " the trackID is : " << track->GetTrackID()<< G4endl;
        GetIT(track)->GetTrackingInfo()->SetLeadingStep(true);
        fLeadingTracks.push_back(track);

    }

    SP->CleanProcessor();
}
//_________________________________________________________________________

void G4ITStepManager::DoIt()

// Call the process having the min step length or just propagate the track on the given time step

// If the track is "leading the step" (ie one of its process has been selected
// as the one having the minimum time step over all tracks and processes),
// it will undergo its selected processes. Otherwise, it will just propagate the track
// on the given time step.

{
    int initialListSize = fpMainList -> size() ;

    for(int i = 0 ; i < initialListSize; i++)
    {
        G4Track* track = fpMainList->back();

        fpMainList->pop_back();
        if(!track)
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "No track was pop back the main track list.";
            G4Exception("G4ITStepManager::DoIt","ITStepManager009",
                        FatalErrorInArgument,exceptionDescription);
        }

        fpStepProcessor -> Stepping(track, fTimeStep) ;

        ExtractDoItData(fpStepProcessor) ;
    }
}
//_________________________________________________________________________

void G4ITStepManager::ExtractDoItData(G4ITStepProcessor* SP)
{

    G4Track* track = SP -> GetTrack() ;
    if(!track)
    {
        SP->CleanProcessor();
        return ;
    }

    G4TrackStatus status  = track->GetTrackStatus() ;

    switch (status)
    {
    case fAlive:
    case fStopButAlive:
    case fSuspend :
    case fPostponeToNextEvent :
    default :
        PushSecondaries(SP);
        fpMainList -> push_front(track);

        break;

    case fStopAndKill :
        PushSecondaries(SP);
        EndTracking(track);
        break;

    case fKillTrackAndSecondaries :
        G4TrackVector* secondaries = SP->GetSecondaries();
        if( secondaries )
        {
            for(size_t i=0; i<secondaries->size(); i++)
            {
                delete (*secondaries)[i];
            }
            secondaries->clear();
        }
        EndTracking(track);
        break;
    }
    SP->CleanProcessor();
}
//_________________________________________________________________________

void G4ITStepManager::PushSecondaries(G4ITStepProcessor* SP)
{
    G4TrackVector* secondaries = SP->GetSecondaries() ;
    if( !secondaries || secondaries -> empty())
    {
        // DEBUG
        //      G4cout << "NO SECONDARIES !!! " << G4endl;
        return ;
    }

    // DEBUG
    //    G4cout << "There are secondaries : "<< secondaries -> size() << G4endl ;

    G4TrackVector::iterator secondaries_i  = secondaries->begin() ;

    for( ; secondaries_i != secondaries->end() ; secondaries_i++)
    {
        G4Track* secondary = *secondaries_i ;
        _PushTrack(secondary);
    }
}

//_________________________________________________________________________

void G4ITStepManager::ComputeTrackReaction()
{
    if(fReactingTracks.empty())
    {
        return ;
    }

    if(fInteractionStep == false)
    {
        fpModelProcessor->FindReaction(&fReactingTracks,
                                       fTimeStep,
                                       fPreviousStepTime,
                                       fReachedUserTimeLimit);
        // TODO
        // A ne faire uniquement si le temps choisis est celui calculé par le time stepper
        // Sinon utiliser quelque chose comme : fModelProcessor->FindReaction(&fMainList);

        std::vector<G4ITReactionChange*> * reactionInfo_v = fpModelProcessor -> GetReactionInfo();
        std::vector<G4ITReactionChange*> :: iterator reactionInfo_i = reactionInfo_v->begin();

        for( ; reactionInfo_i != reactionInfo_v->end() ; reactionInfo_i++)
        {
            G4ITReactionChange* changes = (*reactionInfo_i);
            G4Track* trackA = const_cast<G4Track*>(changes->GetTrackA());
            G4Track* trackB = const_cast<G4Track*>(changes->GetTrackB());

            if(trackA == 0 || trackB == 0 || trackA->GetTrackStatus() == fStopAndKill ||
                    trackB->GetTrackStatus() == fStopAndKill ) continue ;

            G4int nbSecondaries = changes->GetNumberOfSecondaries();

            if(fpUserReactionAction)
            {
                const G4TrackFastVector* productsVector = changes->GetfSecondary();
                fpUserReactionAction->UserReactionAction(*trackA, *trackB, *productsVector, nbSecondaries);
            }

#ifdef G4VERBOSE
            if(fVerbose)
            {
                G4cout << "At time : " << setw(7) << G4BestUnit(fGlobalTime,"Time")
                       << " Reaction : "
                       << GetIT(trackA)->GetName() << " (" << trackA->GetTrackID()
                       << ") + "
                       << GetIT(trackB)->GetName() << " (" << trackB->GetTrackID()
                       << ") -> " ;
            }
#endif

            if(nbSecondaries > 0)
            {
                for(int i = 0 ; i < nbSecondaries ; i ++)
                {
#ifdef G4VERBOSE
                    if(fVerbose && i!=0) G4cout << " + " ;
#endif

                    G4Track* secondary = changes->GetSecondary(i);
                    _PushTrack(secondary);
                    GetIT(secondary)->SetParentID(trackA->GetTrackID(),trackB->GetTrackID());

                    if(secondary->GetGlobalTime() - fGlobalTime > fTimeTolerance)
                    {
                        G4ExceptionDescription exceptionDescription ;
                        exceptionDescription
                                << "The time of the secondary should not be bigger than the current global time."
                                << " This may cause synchronization problem. If the process you are using required "
                                << "such feature please contact the developpers."
                                << G4endl
                                << "The global time in the step manager : " << G4BestUnit(fGlobalTime,"Time")
                                << G4endl
                                << "The global time of the track : " << G4BestUnit(secondary->GetGlobalTime(),"Time")
                                << G4endl;

                        G4Exception("G4ITStepManager::ComputeInteractionBetweenTracks","ITStepManager010",
                                    FatalErrorInArgument,exceptionDescription);
                    }

#ifdef G4VERBOSE
                    if(fVerbose)
                        G4cout << GetIT(secondary)->GetName() << " (" << secondary->GetTrackID() <<")";
#endif
                }
            }
            else
            {
#ifdef G4VERBOSE
                if(fVerbose)
                    G4cout<<"No product vector" ;
#endif
            }
#ifdef G4VERBOSE
            if(fVerbose)
                G4cout << G4endl;
#endif
            if(trackA->GetTrackID() == 0 || trackB->GetTrackID() == 0)
            {
                G4Track* track = 0;
                if(trackA->GetTrackID() == 0) track = trackA;
                else track = trackB;

                G4ExceptionDescription exceptionDescription ;
                exceptionDescription << "The problem was found for the reaction between tracks :"
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
                                     <<   trackA->GetParentID() << "\n";
                exceptionDescription << "Parent ID of trackB : "
                                     <<   trackB->GetParentID() << "\n";

                exceptionDescription << "The ID of one of the reaction track was not setup.";
                G4Exception("G4ITStepManager::ComputeInteractionBetweenTracks","ITStepManager011",
                            FatalErrorInArgument,exceptionDescription);
            }

            if(changes->WereParentsKilled())
            {
                // DEBUG
                //                G4cout << "Erasing tracks : "
                //                 << "trackA at time : " << G4BestUnit(trackA->GetGlobalTime() , "Time")
                //                 << "\t trackB at time : "<< G4BestUnit(trackB->GetGlobalTime(), "Time")
                //                 << "\t GlobalTime : " << G4BestUnit(fGlobalTime, "Time")
                //                 << G4endl;

                trackA->SetTrackStatus(fStopAndKill);
                trackB->SetTrackStatus(fStopAndKill);
                if(fpMainList->Holds(trackA))  fpMainList->pop(trackA);
                else                           fSecondaries.pop(trackA);

                if(fpMainList->Holds(trackB))  fpMainList->pop(trackB);
                else                           fSecondaries.pop(trackB);
                EndTracking(trackA);
                EndTracking(trackB);
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

    fReactingTracks.clear();
}
//_________________________________________________________________________

void G4ITStepManager::MergeSecondariesWithMainList()
{
    fSecondaries.transferTo(fpMainList);
}

//_________________________________________________________________________

void G4ITStepManager::AddTrackID(G4Track* track)
{
    if(fNbTracks == 0) fNbTracks = -1;
    track->SetTrackID(fNbTracks);
    fNbTracks--;
}

//_________________________________________________________________________

void G4ITStepManager::PushTrack(G4Track* track)
{
    if(fRunning)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription
                << "G4ITStepManager::PushTrack : You are trying to push tracks while the ITStepManager is running";
        G4Exception("G4ITStepManager::PushTrack","ITStepManager012",
                    FatalErrorInArgument,exceptionDescription);
    }
    _PushTrack(track);
}
//_________________________________________________________________________

void G4ITStepManager::_PushTrack(G4Track* track)
{
    G4double globalTime = track->GetGlobalTime() ;

    if(track->GetTrackID() == 0)
    {
        // Set track ID
        AddTrackID(track);
    }

#ifdef G4VERBOSE
    if(fVerbose > 2)
    {
        G4cout << "Pushing a track : trackID " << track->GetTrackID() << G4endl;
        G4cout << "The time of the ITStepManager: " << G4BestUnit(fGlobalTime,"Time") << G4endl;
        G4cout << "The track's time: " << G4BestUnit(track->GetGlobalTime(),"Time") << G4endl;
    }
#endif

    if(!fRunning)
    {
        if(globalTime < fGlobalTime)
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription
                    << "You are trying to push a track with a global time"
                    << " inferior to the current simulation time." << G4endl
                    << "The time is going back : " << G4endl
                    << "The time in the step manager : " << G4BestUnit(fGlobalTime,"Time")<< G4endl
                    << "The time of the track : " << G4BestUnit(globalTime,"Time") << G4endl
                    << "(ITStepManager is not yet running)" << G4endl;

            G4Exception("G4ITStepManager::_PushTrack","ITStepManager014",
                        FatalErrorInArgument,exceptionDescription);
        }

        // Push the track to the rigth track list :
        // If the track time is the same as the main track list,
        // it will be push to the main track list
        // otherwise, it will be pushed to the delayed track list.
        if(fpMainList && fpMainList->empty())
        {
            PushDelayed(track, globalTime);
        }
        else
        {
            if(fpMainList && globalTime == fGlobalTime)
            {
                fpMainList->push_back(track);
            }
            else
            {
                PushDelayed(track, globalTime);
            }
        }
    }
    else
    {
        double timeDifference = globalTime - fGlobalTime ;

        if(timeDifference < -1*fTimeTolerance)
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription
                    << "You are trying to push a track with a global time"
                    << " inferior to the current simulation time." << G4endl
                    << "The time is going back : " << G4endl
                    << "The time in the step manager : " << G4BestUnit(fGlobalTime,"Time")<< G4endl
                    << "The time of the track : " << G4BestUnit(globalTime,"Time") << G4endl
                    << "(ITStepManager is running)" << G4endl;

            G4Exception("G4ITStepManager::_PushTrack","ITStepManager015",
                        FatalErrorInArgument,exceptionDescription);
        }

        // Push the track to the rigth track list :
        // If the track time is the same as the main track list,
        // it will be push to the secondary list
        // otherwise, it will be pushed to the delayed track list.
        if(fabs(timeDifference) < fTimeTolerance)
        {
            fSecondaries.push_back(track);
        }
        else // globalTime < fGlobalTime already taken into account above
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription
                    << "While running you cannot push a track"
                    << " with a bigger global time than the current global time" << G4endl
                    << "The time in the step manager : " << G4BestUnit(fGlobalTime,"Time")<< G4endl
                    << "The time of the track : " << G4BestUnit(globalTime,"Time") << G4endl
                    << "(ITStepManager is running)" << G4endl;

            G4Exception("G4ITStepManager::_PushTrack","ITStepManager016",
                        FatalErrorInArgument,exceptionDescription);
            // PushDelayed(track, globalTime);
        }
    }
}
//_________________________________________________________________________

void G4ITStepManager::PushDelayed(G4Track* track, const G4double& globalTime)
{
#ifdef G4VERBOSE
    if(fVerbose > 2)
    {
        G4cout << "Pushing a delayed track : trackID " << track->GetTrackID() << G4endl;
        G4cout << "The time of the ITStepManager: " << G4BestUnit(fGlobalTime,"Time") << G4endl;
        G4cout << "The track's time: " << G4BestUnit(track->GetGlobalTime(),"Time") << G4endl;
    }
#endif

    std::map<double,G4TrackList* > :: iterator
            fDelayedList_i = fDelayedList.find(globalTime) ;

    if(fDelayedList_i == fDelayedList.end())
    {

        G4TrackList* newList = new G4TrackList ;
        newList -> push_back(track);
        fDelayedList[globalTime] = newList ;
    }
    else
    {
        fDelayedList_i->second-> push_back(track);
    }
}
//_________________________________________________________________________

void G4ITStepManager::KillTracks()
{
    fToBeKilledList.erase(fToBeKilledList.begin(), fToBeKilledList.end());
}
//_________________________________________________________________________

void G4ITStepManager::EndTracking(G4Track* trackToBeKilled)
{
    fpTrackingManager->EndTracking(trackToBeKilled);
    fToBeKilledList.push_back(trackToBeKilled);
}
//_________________________________________________________________________

void G4ITStepManager::EndTracking()
{
    if(fpMainList && !fpMainList->empty())
    {
        G4TrackList::iterator      fpMainList_i = fpMainList->begin();
        for( ; fpMainList_i != fpMainList->end() ; fpMainList_i++)
        {
            fpTrackingManager->EndTracking(*fpMainList_i);
        }
    }

    if(!fSecondaries.empty())
    {
        G4TrackList::iterator      fSecondaries_i = fSecondaries.begin();
        for( ; fSecondaries_i != fSecondaries.end() ; fSecondaries_i++)
        {
            fpTrackingManager->EndTracking(*fSecondaries_i);
        }
    }
    if(fpWaitingList && !fpWaitingList->empty())
    {
        G4TrackList::iterator      fpWaitingList_i = fpWaitingList->begin();
        for( ; fpWaitingList_i != fpWaitingList->end() ; fpWaitingList_i++)
        {
            fpTrackingManager->EndTracking(*fpWaitingList_i);
        }
    }
}

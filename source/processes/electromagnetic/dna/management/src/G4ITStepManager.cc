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
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
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
#include <stdlib.h>
#include "G4ITTrackingManager.hh"
#include "G4TrackingInformation.hh"

using namespace std;

G4ITStepManager* G4ITStepManager::fgStepManager = 0 ;

template<typename T>
inline bool IsInf(T value)
{
    return std::numeric_limits<T>::has_infinity &&
            value == std::numeric_limits<T>::infinity();
}
//_________________________________________________________________________

G4ITStepManager* G4ITStepManager::Instance()
{
    if(!fgStepManager) fgStepManager = new G4ITStepManager();
    return fgStepManager ;
} 
//_________________________________________________________________________

G4ITStepManager::G4ITStepManager()
{
    //    fNbProc = get_nprocs ();
    fNbProc = 1 ;
    // DEBUG
    //    G4cout <<"Number of processors detected : " << fNbProc << G4endl;

    fpMainList = 0 ;
    fpWaitingList = 0 ;

    fpUserTimeSteps = 0;

    fMinTimeStep = DBL_MAX ;
    fComputeTimeStep = false;
    fComputeReaction = false;
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

    fpUserITAction = 0;
    fpTrackingManager = 0;

    fNbTracks = -1;
    fpModelHandler = new G4ITModelHandler();
    fpTrackingManager = new G4ITTrackingManager();
}
//_________________________________________________________________________

void G4ITStepManager::ClearList()
{
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

    std::map<double,G4TrackList* >::iterator fDelayedList_i = fDelayedList.begin() ;

    for(; fDelayedList_i != fDelayedList.end() ; fDelayedList_i++)
    {
        if(fDelayedList_i->second)
            delete (fDelayedList_i->second);
        fDelayedList_i->second = 0 ;
    }
    fDelayedList.clear();
    fNbTracks = -1;
}
//_________________________________________________________________________

G4ITStepManager::~G4ITStepManager()
{
    delete fpMasterStepProcessor ;
    delete fpMasterModelProcessor ;
    delete G4ITTypeManager::Instance();
    G4ITStepManager::Instance() -> ClearList();
    delete fpTrackingManager;
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

    fpMasterModelProcessor = new G4ITModelProcessor(fNbProc-1);
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
    fMinTimeStep = DBL_MAX ;
    fGlobalTime = -1 ;
    fInteractionStep = true ;

    fNbSteps = 0;
}
//_________________________________________________________________________

void G4ITStepManager::Process()
{
    fpTrackingManager->Initialize();
    fpMasterModelProcessor->Initialize();
    fpMasterStepProcessor->Initialize();

    // ___________________
    fRunning = true ;
    if( ! fDelayedList.empty())     SynchronizeTracks() ;
    DoProcess() ;

    EndTracking();

    // ___________________
    fRunning = false;
}
//_________________________________________________________________________

void G4ITStepManager::SynchronizeTracks()
{
    if(fpWaitingList != 0)
    {
        __Exception_Origin__
        G4String exceptionCode ("ITStepManager002");
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "There is a waiting track list (fpWaitingList != 0).";
        exceptionDescription << " When G4ITStepManager::SynchronizeTracks() is called, ";
        exceptionDescription << "no more tracks should remain in the fpWaitingList.";
        G4Exception(exceptionOrigin.data(),exceptionCode.data(),
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
    
    fpWaitingList    = fpMainList ;
    fTmpGlobalTime  = fGlobalTime ;
    fTmpEndTime     = fEndTime ;

    std::map<double,G4TrackList* >::iterator fDelayedList_i = fDelayedList.begin();
    std::map<double,G4TrackList* >::iterator nextfDelayedList_i;
    if(fDelayedList_i != fDelayedList.end())
    {
        // Retrieve the first delayed list
        fpMainList       = fDelayedList_i->second ;
        fGlobalTime     = fDelayedList_i->first;

        do
        {
            if((fDelayedList_i->second) == 0)   continue ;
            if(fDelayedList_i->second->empty()) continue ;

            if(fDelayedList_i->first > fTmpGlobalTime)
            {

                G4cerr << "Next time to consider = " << fDelayedList_i->first / picosecond << G4endl;
                G4cerr << "G4ITStepManager::fTmpGlobalTime (should big the biggest time in the simulation) = "
                       << fTmpGlobalTime / picosecond<< G4endl;

                __Exception_Origin__
                G4String exceptionCode ("ITStepManager003");
                G4ExceptionDescription exceptionDescription ;
                exceptionDescription << "The next time in delayed tracks is bigger than the maximum allowed simulation time.";
                G4Exception(exceptionOrigin.data(),exceptionCode.data(),
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
    G4cout << "G4ITStepManager::DoProcess()"<<G4endl;
    while(fGlobalTime <= fEndTime &&  ! (fpMainList -> empty()) )
    {
        Stepping() ;
    }

    G4cout << "!!!******** Ending at time : " << fGlobalTime / microsecond << "[microsecond]"<< G4endl;
}
//_________________________________________________________________________

void G4ITStepManager::Stepping()
{

    // DEBUG
    //G4cout << "In G4ITStepManager::Stepping()" << G4endl;

    // Before anything move PostStepPoint to PreStepPoint for every tracks

    fMinTimeStep = DBL_MAX ;
    fInteractionStep = false ;
    fReachedUserTimeLimit = false;

    FindUserPreDefinedTimeStep();

    if(fComputeTimeStep)
    {
        CalculateMinStep() ; // => at least N (N = nb of tracks) loops
    }
    else
    {
        fMinTimeStep = fDefinedMinTimeStep ;
    }

    // DEBUG
//        G4cout << "********* HERE ******** " << G4endl;
//        G4cout << "\33[1;31m"
//                 << "CalculateMinStep :"
//                 << " !!!! ***** fMinTimeStep : "
//                 << fMinTimeStep /picosecond
//                 << " [ps]"
//                 << "\33[0m"
//                 << G4endl;
//        G4cout << "ComputeInteractionLength()" << G4endl;


    if(fMinTimeStep > 0)
    {

        ComputeInteractionLength(); // => at least N loops
        // All process returns the physical step of interaction
        // using the approximation : E = cste along the step
        // Only the transportation calculates the corresponding
        // time step
        // DEBUG
//        G4cout << "(ComputeInteractionLength) :"  << " !!!! **** Time Step : "     << fMinTimeStep /picosecond << " [ps]"<< G4endl;
//        G4cout << " fInteractionStep : " << fInteractionStep << G4endl;


        fReachedUserTimeLimit  = (fMinTimeStep == fDefinedMinTimeStep)    ? true : false;
        // DEBUG
//        G4cout << " !!!! **** Global time : "   << fGlobalTime  /picosecond << " [ps]"<< G4endl;
//        G4cout << "DoIt()" << G4endl;

        DoIt(); // should redo transportation ASIL

        fGlobalTime += fMinTimeStep ;

        MergeSecondariesWithMainList();
        KillTracks();
    }
    else
    {
        fReachedUserTimeLimit = (fMinTimeStep == fDefinedMinTimeStep)    ?  true : false;
    }
    // DEBUG
    //    G4cout << "ComputeInteractionBetweenTracks()" << G4endl;

    ComputeInteractionBetweenTracks();

    MergeSecondariesWithMainList();

    KillTracks();

    fNbSteps++;

    // DEBUG
    //    G4cout  << "___________________________________" << G4endl;
    //    G4cout  << "*** End Of Step ***" << G4endl;
    //    G4cout  << "Previous global time : " << (fGlobalTime- fMinTimeStep)
    //                                    / picosecond
    //            << "[ps]\t Chosen time :"
    //            << fMinTimeStep / picosecond<<G4endl;
    //    G4cout  << "Next global time : " << fGlobalTime
    //               / picosecond << G4endl;
    //    G4cout  << "___________________________________" << G4endl;
    

}
//_________________________________________________________________________

void G4ITStepManager::FindUserPreDefinedTimeStep()
{
    if(fUsePreDefinedTimeSteps)
    {
        if(! fpUserTimeSteps)
        {
            __Exception_Origin__
            G4String exceptionCode ("ITStepManager004");
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "You are asking to use user defined steps but you did not give any.";
            G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                        FatalErrorInArgument,exceptionDescription);
        }
        map<double, double>::iterator fpUserTimeSteps_i   = fpUserTimeSteps->upper_bound(fGlobalTime) ;
        map<double, double>::iterator fpUserTimeSteps_low = fpUserTimeSteps->lower_bound(fGlobalTime) ;

        // DEBUG
        //        G4cout << "fGlobalTime : " << fGlobalTime /picosecond  << G4endl;
        //        G4cout << "fpUserTimeSteps_i : "
        //        <<"<"<<fpUserTimeSteps_i->first / picosecond<<", "<< fpUserTimeSteps_i->second / picosecond<<">"
        //        << "\t fpUserTimeSteps_low : "
        //        <<"<"<<fpUserTimeSteps_low->first / picosecond<<", "<< fpUserTimeSteps_low->second / picosecond<<">"
        //        << G4endl;


        if(fpUserTimeSteps_i == fpUserTimeSteps->end())
        {
            // DEBUG
            //            G4cout << "1" << G4endl;
            fpUserTimeSteps_i --;
        }
        else if(fabs(fGlobalTime - fpUserTimeSteps_low->first) < fGlobalTime*1/100. )
        {
            // Case : fGlobalTime = X picosecond
            // and fpUserTimeSteps_low->first = X picosecond
            // but the precision is not good enough
            // DEBUG
            //            G4cout << "2" << G4endl;
            fpUserTimeSteps_i = fpUserTimeSteps_low;
        }
        else if(fpUserTimeSteps_i == fpUserTimeSteps_low)
        {
            // "Normal" cases
            // DEBUG
            //            G4cout << "3" << G4endl;
            fpUserTimeSteps_i--;
        }
        else
        {
            // DEBUG
            //            G4cout << "4" << G4endl;
            fpUserTimeSteps_i = fpUserTimeSteps_low;
            //fDefinedMinTimeStep = fpUserTimeSteps_low->second ;
        }

        fDefinedMinTimeStep = fpUserTimeSteps_i->second ;

        // DEBUG
        //        G4cout << "fGlobalTime : " << fGlobalTime /picosecond << "\t chosen step time : " << fDefinedMinTimeStep /picosecond<< G4endl;
        //        G4cout << "fpUserTimeSteps_i : "<< fpUserTimeSteps_i->second / picosecond << "\t fpUserTimeSteps_low : "<< fpUserTimeSteps_low->second / picosecond<< G4endl;
    }
    else fDefinedMinTimeStep = 0 ;
}
//_________________________________________________________________________

void G4ITStepManager::CalculateMinStep()
{

    if(fpMasterModelProcessor == 0)
    {
        __Exception_Origin__
        G4String exceptionCode ("ITStepManager005");
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "There is no G4ITModelProcessor to hande IT reaction. ";
        exceptionDescription << "You probably did not initialize the G4ITStepManager. ";
        exceptionDescription << "Just do G4ITStepManager::Instance()->Initialize(); ";
        exceptionDescription << " but only after initializing the run manager.";
        G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                    FatalErrorInArgument,exceptionDescription);
    }

    fpMasterModelProcessor -> InitializeStepper(fGlobalTime, fDefinedMinTimeStep) ;

    G4TrackList::iterator      fpMainList_i = fpMainList->begin();
    for( ; fpMainList_i != fpMainList->end() ; fpMainList_i++)
    {
        G4Track * track  = *fpMainList_i ;

        if(track == 0)
        {
            __Exception_Origin__
            G4String exceptionCode ("ITStepManager006");
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "No track found.";
            G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                        FatalErrorInArgument,exceptionDescription);
        }

        if(track->GetPosition() == G4ThreeVector())
        {
            G4cout << "track at origin" << G4endl;
            G4cout << "track Id : " << track->GetTrackID() << G4endl;
        }

        G4TrackStatus trackStatus = track->GetTrackStatus();
        if(trackStatus == fStopAndKill || trackStatus == fStopButAlive)
        {
            continue ;
        }

        while (fpModelProcessor->Occupied())
        {
            fpModelProcessor = fpModelProcessor -> GetNext() ;
        };

        ExtractTimeStepperData(fpModelProcessor) ;

        //DEBUG
        //        if(fUsePreDefinedTimeSteps)
        //        {
        //                G4cout << "Calculate step for track : "
        //                 << track->GetTrackID()
        //                 << " in stepProcessor : " << fpStepProcessor
        //                 << G4endl;

        // DEBUG
        //  G4cout << "send to model processor" << G4endl;
        fpModelProcessor->CalculateStep(track, fDefinedMinTimeStep);

        //        } //DEBUG

    }

    // DEBUG
    // G4cout << "done" << G4endl;

    Join(fpModelProcessor);


    ExtractRemainingData(&G4ITStepManager::ExtractTimeStepperData, fpModelProcessor) ;

    // DEBUG
    //        G4cout << "done 2 | time : " << fGlobalTime /picosecond << G4endl;
    //     std::map<G4Track*, std::vector<G4Track*>* > ::iterator fReactingTracks_i ;
    //
    //     for(fReactingTracks_i = fReactingTracks.begin() ; fReactingTracks_i != fReactingTracks.end()
    //             ; fReactingTracks_i ++)
    //     {
    //         G4Track* trackA = fReactingTracks_i->first;
    //         std::vector<G4Track*>* trackB_vector = fReactingTracks_i->second ;
    //         std::vector<G4Track*>::iterator trackB_i = trackB_vector->begin();
    //         G4Track* trackB = 0 ;
    //
    //         G4cout << "track A : "  << trackA->GetTrackID() << G4endl;
    //
    //         for(; trackB_i != trackB_vector->end() ; trackB_i++)
    //         {
    //             trackB = *trackB_i;
    //             G4cout << "\t\t track B : "  << trackB->GetTrackID() << G4endl;
    //
    //         }
    //
    //     }

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
    // DEBUG
    //    else
    //    {
    //      G4cout << "track address : " << track << G4endl;
    //      G4cout << "Particle name : " << track->GetParticleDefinition()->GetParticleName() << G4endl;
    //    }

    const std::vector<std::vector<G4VITModel*> >* model = MP->GetCurrentModel();

    for(int i = 0 ; i < (int) model->size() ; i++)
    {
        for(int j = 0 ; j < (int) (*model)[i].size() ; j++)
        {
            G4VITModel* mod = (*model)[i][j];

            if(mod == 0)
            {
                continue;
            }

            G4VITTimeStepper* stepper (mod -> GetTimeStepper()) ;

            G4double sampledMinTimeStep (stepper -> GetSampledMinTimeStep()) ;
            G4TrackVectorHandle reactants (stepper -> GetReactants());

            if(sampledMinTimeStep < fMinTimeStep)
            {
                if(!reactants)
                {
                    // FIXME : avoid the line below
                    if(track->GetParticleDefinition()->GetParticleName() == "H2O") continue;
                    G4cerr<< "reactants == 0"<< G4endl;
                    G4cerr<< "sampledMinTimeStep : " << sampledMinTimeStep<< G4endl;

                    __Exception_Origin__
                    G4String exceptionCode ("ITStepManager007");
                    G4ExceptionDescription exceptionDescription ;
                    exceptionDescription << "No reactants were registered.";
                    G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                                FatalErrorInArgument,exceptionDescription);

                    continue ;
                }

                fMinTimeStep = sampledMinTimeStep ;
                fReactingTracks.clear();
                fReactingTracks.insert(make_pair(track, reactants));

                stepper -> ResetReactants();

                // DEBUG
                //                G4cout << "Reactants ptr : " << stepper -> GetReactants() << G4endl;
                //                for(G4int i = 0 ; i < fReactingTracks[track]->size() ; i++)
                //                {
                //                	G4cout << "\t\t track ID in fReactingTracks : "
                //                	<< (*fReactingTracks[track])[i]->GetTrackID() << G4endl;
                //
                //                }
                //
                //                 G4cout << "-------------------------" << G4endl;
                //                 G4cout << "Stepper has found reaction for "
                //                 << GetMolecule(track)->GetName()
                //                 << " (" <<  track->GetTrackID()<<") | at time "
                //                 << (fGlobalTime + sampledMinTimeStep) / picosecond
                //                 << " corresponding TS : " << sampledMinTimeStep /picosecond << "[ps]"
                //                 << G4endl;
                //                 for(G4int i = 0 ; i < (G4int) fReactingTracks[track]->size() ; i++)
                //                 {
                //                     G4cout <<"\t \t \t"<<  GetMolecule((*fReactingTracks[track])[i])->GetName()
                //                     << " ID : "<< (*fReactingTracks[track])[i]->GetTrackID()
                //                     << " distance : "
                //                     << (track->GetPosition() - (*fReactingTracks[track])[i]->GetPosition()).mag()
                //                     << G4endl;
                //
                //                     if(track == (*fReactingTracks[track])[i])
                //                     {
                //                         G4cout << "BIG PB" << G4endl;
                //                         G4Exception();
                //                     }
                //                 }
                //
            }
            else if(fMinTimeStep == sampledMinTimeStep )
            {

                if(!reactants)
                {
                    G4cerr<< "reactants == 0"<< G4endl;
                    G4cerr<< "sampledMinTimeStep : " << sampledMinTimeStep<< G4endl;

                    __Exception_Origin__
                    G4String exceptionCode ("ITStepManager008");
                    G4ExceptionDescription exceptionDescription ;
                    exceptionDescription << "No reactants were registered.";
                    G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                                FatalErrorInArgument,exceptionDescription);
                }

                fReactingTracks.insert(make_pair(track, reactants));

                stepper -> ResetReactants();

                // DEBUG
                //
                //                G4cout << "-------------------------" << G4endl;
                //                 G4cout << "Stepper has found reaction for "
                //                 << GetMolecule(track)->GetName()
                //                 << " (" <<  track->GetTrackID()<<") | at time "
                //                 << (fGlobalTime + sampledMinTimeStep) / picosecond
                //                 << " corresponding TS : " << sampledMinTimeStep /picosecond << "[ps]"
                //                 << G4endl;
                //                 for(G4int i = 0 ; i < (G4int) fReactingTracks[track]->size() ; i++)
                //                 {
                //                     G4cout <<"\t \t \t"<<  GetMolecule((*fReactingTracks[track])[i])->GetName()
                //                     << " ID : "<< (*fReactingTracks[track])[i]->GetTrackID() << G4endl;
                //
                //                     if(track == (*fReactingTracks[track])[i])
                //                     {
                //                         G4cout << "BIG PB" << G4endl;
                //                         G4Exception();
                //                     }
                //                 }
                //
            }
            else
            {
                if(fabs(fMinTimeStep-sampledMinTimeStep) <= 1e-12)
                {
                    G4cerr << "fMinTimeStep : " << fMinTimeStep / picosecond
                           << "[ps]"<<G4endl;
                    G4cerr << "sampledMinTimeStep : " << sampledMinTimeStep/ picosecond
                           << "[ps]"<<G4endl;

                    // TODO
                    // Find a exception description
                    __Exception_Origin__
                    G4String exceptionCode ("ITStepManager009");
                    G4ExceptionDescription exceptionDescription ;
                    exceptionDescription << " ... ";
                    G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                                FatalErrorInArgument,exceptionDescription);
                }
                stepper -> ResetReactants();
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
    else if( SP -> GetInteractionTime() < fMinTimeStep)
    {
        if(fInteractionStep == false)
        {
            fInteractionStep = true ;
        }
        else
        {
            std::vector<G4Track*>::iterator fLeadingTracks_i = fLeadingTracks.begin();

            while(fLeadingTracks_i != fLeadingTracks.end())
            {
                G4Track* track = *fLeadingTracks_i;
                if(!track)
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
//               << SP -> GetInteractionTime() << " against fMinTimeStep : "
//               << fMinTimeStep << " the trackID is : " << track->GetTrackID() << G4endl;

        fMinTimeStep =  SP -> GetInteractionTime() ;

        GetIT(track)->GetTrackingInfo()->SetLeadingStep(true);
        fLeadingTracks.push_back(track);
    }
    else if(fMinTimeStep == SP -> GetInteractionTime() )
    {
        if(! fInteractionStep)
        {
            // Attention fMinTimeStep == fMinReactionTime
            fInteractionStep = true ;
        }

//        G4cout << "Will set leading step to true for time  :"
//               << SP -> GetInteractionTime() << " against fMinTimeStep : "
//               << fMinTimeStep << " the trackID is : " << track->GetTrackID()<< G4endl;
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

            __Exception_Origin__
            G4String exceptionCode ("ITStepManager010");
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "No track was pop back the main track list.";
            G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                        FatalErrorInArgument,exceptionDescription);
        }
        fpStepProcessor -> Stepping(track, fMinTimeStep) ;

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
        secondary->SetTrackID(fNbTracks);
        fNbTracks--;
        fSecondaries.push_back(secondary);
    }

}

//_________________________________________________________________________

void G4ITStepManager::ComputeInteractionBetweenTracks()

{
    if(fReactingTracks.empty())
    {
        return ;
    }

    if(fInteractionStep == false)
    {
        fpModelProcessor->FindReaction(&fReactingTracks,
                                       fMinTimeStep,
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

            if(fpUserITAction)
            {
                const G4ConstTrackFastVector* trackVector = (G4ConstTrackFastVector*) changes->GetfSecondary();
                fpUserITAction->UserReactionAction(*trackA, *trackB, *trackVector, nbSecondaries);
            }

            G4cout << "Reaction : "
                   << GetIT(trackA)->GetName() << " (" << trackA->GetTrackID()
                   << ") + "
                   << GetIT(trackB)->GetName() << " (" << trackB->GetTrackID()
                   << ") -> " ;

            if(nbSecondaries > 0)
            {
                for(int i = 0 ; i < nbSecondaries ; i ++)
                {
                    if(i!=0) G4cout << " + " ;
                    G4Track* secondary = changes->GetSecondary(i);

                    if(secondary->GetTrackID() == 0)
                    {
                        if(fNbTracks == 0) fNbTracks = -1;
                        secondary->SetTrackID(fNbTracks);
                        fNbTracks--;
                    }
                    fpMainList->push_back(secondary);

                    if(secondary->GetGlobalTime() - fGlobalTime > (1-1/100)*fGlobalTime)
                    {
                        G4cerr << secondary->GetGlobalTime()  << G4endl;
                        G4cerr << fGlobalTime  << G4endl;

                        __Exception_Origin__
                        G4String exceptionCode ("ITStepManager011");
                        G4ExceptionDescription exceptionDescription ;
                        exceptionDescription << "The time of the secondary is going back ";
                        exceptionDescription << "(ie the time of the secondary in the reference frame is less than ";
                        exceptionDescription << "the global time of the simulation).";
                        G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                                    FatalErrorInArgument,exceptionDescription);
                    }

                    G4cout << GetIT(secondary)->GetName() << " (" << secondary->GetTrackID() <<")";
                }
            }
            else
            {
                G4cout<<"No product vector" ;
            }

            G4cout << " | At time : " << fGlobalTime /picosecond << G4endl;

            if(trackA->GetTrackID() == 0 || trackB->GetTrackID() == 0)
            {
                G4cerr<<G4endl;
                G4cerr << "***** PROBLEM *****" << G4endl;
                G4Track* track = 0;
                if(trackA->GetTrackID() == 0) track = trackA;
                else track = trackB;

                G4cerr << "The problem was found for the reaction between tracks :"
                       << trackA->GetParticleDefinition()->GetParticleName() << " ("
                       << trackA->GetTrackID() << ") & "
                       << trackB->GetParticleDefinition()->GetParticleName() << " ("
                       << trackB->GetTrackID() << ")."<< G4endl;

                if(track->GetStep() == 0)
                {
                    G4cerr << "Also no step was found"
                           << " ie track->GetStep() == 0 "
                           << G4endl;
                }

                G4cerr << "Parent ID of trackA : "
                       <<   trackA->GetParentID() << G4endl;
                G4cerr << "Parent ID of trackB : "
                       <<   trackB->GetParentID() << G4endl;


                __Exception_Origin__
                G4String exceptionCode ("ITStepManager012");
                G4ExceptionDescription exceptionDescription ;
                exceptionDescription << "The ID of one of the reaction track was not setup.";
                G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                            FatalErrorInArgument,exceptionDescription);
            }

            if(changes->WereParentsKilled())
            {
                // DEBUG
                //                G4cout << "Erasing tracks : "
                //                 << "trackA at time : " << trackA->GetGlobalTime() / picosecond
                //                 << "\t trackB at time : "<< trackB->GetGlobalTime() / picosecond
                //                 << "\t GlobalTime : " << fGlobalTime / picosecond
                //                 << G4endl;

                trackA->SetTrackStatus(fStopAndKill);
                trackB->SetTrackStatus(fStopAndKill);
                fpMainList->pop(trackA);
                fpMainList->pop(trackB);
                EndTracking(trackA);
                EndTracking(trackB);
            }

            delete changes;
        }

        fReactingTracks.erase(fReactingTracks.begin(), fReactingTracks.end());
        reactionInfo_v->clear();
    }
    // DEBUG
    //else
    //{
    //    G4cout << "fInteractionStep == true" << G4endl ;
    //}
}
//_________________________________________________________________________

//void G4ITStepManager::CallPreTrackingProcessesForSecondaries()
//{
//    fpMasterModelProcessor -> InitializePreTrackingProc() ;

//    G4TrackList::iterator fSecondaries_i = fSecondaries.begin();

//    while(fSecondaries_i!= fSecondaries.end())
//    {
//        ExtractPreTrackingProcData(fpModelProcessor) ;

//        G4Track* track = fSecondaries.back();
//        // Make Action
//        fSecondaries_i ++ ;
//    }
//}
//_________________________________________________________________________

//void G4ITStepManager::ExtractPreTrackingProcData(G4ITModelProcessor* SP)
//{
//    G4Track* track = SP -> GetTrack() ;
//    if(!track) return ;

//    SP->CleanProcessor();
//}
//_________________________________________________________________________

void G4ITStepManager::MergeSecondariesWithMainList()
{
    fSecondaries.transferTo(fpMainList);
}
//_________________________________________________________________________

void G4ITStepManager::PushTrack(G4Track* track)
{
    if(fRunning)
    {
        __Exception_Origin__
        G4String exceptionCode ("ITStepManager013");
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription
        << "G4ITStepManager::PushTrack : You are trying to push tracks while the ITStepManager is running";
        G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                    FatalErrorInArgument,exceptionDescription);
    }
    _PushTrack(track);
}
//_________________________________________________________________________

void G4ITStepManager::_PushTrack(G4Track* track)
{
    G4double globalTime = track->GetGlobalTime() ;

    if(globalTime < fGlobalTime)
    {
        G4cerr<< "Erreur : G4ITStepManager le temps recul"<<G4endl;
        G4cerr<<"fGlobalTime : " << fGlobalTime /picosecond << "[ps]" << G4endl;
        G4cerr<<"globalTime : " << globalTime   /picosecond << "[ps]" << G4endl;

        __Exception_Origin__
        G4String exceptionCode ("ITStepManager014");
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription
        << "G4ITStepManager::PushTrack : You are trying to push tracks with a global time inferior to the current simulation time ";
        G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                    FatalErrorInArgument,exceptionDescription);
    }

    // Set track ID
    if(track->GetTrackID() == 0)
    {
        if(fNbTracks == 0) fNbTracks = -1;
        track->SetTrackID(fNbTracks);
        fNbTracks--;
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
//_________________________________________________________________________

void G4ITStepManager::PushDelayed(G4Track* track, const G4double& globalTime)
{

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

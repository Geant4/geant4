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

#include "G4ITStepProcessor.hh"
#include "G4UImanager.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4TransportationManager.hh"
// #include "G4VSensitiveDetector.hh"    // Include from 'hits/digi'
#include "G4GeometryTolerance.hh"
#include "G4ParticleTable.hh"
#include "G4ITTrackingManager.hh"
#include "G4TrackingInformation.hh"
#include "G4IT.hh"

#include "G4VITProcess.hh"
#include "G4VProcess.hh"
#include "G4ITTransportation.hh"

#include <iomanip>              // Include from 'system'
#include <vector>               // Include from 'system'

using namespace std;

static const size_t SizeOfSelectedDoItVector=100;
//____________________________________________________________________________________

G4ITStepProcessor::G4ITStepProcessor()
{
    verboseLevel = 0 ;
//    fpUserSteppingAction = 0 ;
    fpTrackingManager = 0;
    fpNavigator = 0;
    fpTouchableHandle = 0;
    CleanProcessor();
}
//____________________________________________________________________________________

void G4ITStepProcessor::Initialize(void*)
{
    CleanProcessor();
    ActiveOnlyITProcess();

    SetNavigator(G4TransportationManager::GetTransportationManager()
                 ->GetNavigatorForTracking());

    physIntLength = DBL_MAX;
    kCarTolerance = 0.5*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}
//______________________________________________________________________________

G4ITStepProcessor::~G4ITStepProcessor()
{
    if(fpStep)
    {
        fpStep->DeleteSecondaryVector();
        delete fpStep;
    }

    if(fpSecondary)                      delete fpSecondary;
    if(fpSelectedAtRestDoItVector)       delete fpSelectedAtRestDoItVector;
    if(fpSelectedAlongStepDoItVector)    delete fpSelectedAlongStepDoItVector;
    if(fpSelectedPostStepDoItVector)     delete fpSelectedPostStepDoItVector;
//    if(fpUserSteppingAction)             delete fpUserSteppingAction;
}
//______________________________________________________________________________
// should not be used
G4ITStepProcessor::G4ITStepProcessor(const G4ITStepProcessor& rhs)
{
    verboseLevel = rhs.verboseLevel ;
//    fpUserSteppingAction = 0 ;
    fpTrackingManager = 0;
    fpNavigator = 0;
    fpTouchableHandle = 0;
    CleanProcessor();
}
//______________________________________________________________________________

G4ITStepProcessor& G4ITStepProcessor::operator=(const G4ITStepProcessor& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}
// ******************************************************************

void G4ITStepProcessor::ActiveOnlyITProcess()
{
#ifdef debug
    G4cout<<"G4ITStepProcessor::CloneProcesses: is called"<<G4endl;
#endif

    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable->GetIterator();

    theParticleIterator->reset();
    // TODO : Ne faire la boucle que sur les IT **** !!!
    while( (*theParticleIterator)() )
    {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pm= particle->GetProcessManager();

        if(!pm)
        {
            G4cerr << "ERROR - G4ITStepProcessor::GetProcessNumber()" << G4endl
                   << "        ProcessManager is NULL for particle = "
                   << particle->GetParticleName() << ", PDG_code = "
                   << particle->GetPDGEncoding() << G4endl;
            G4Exception("G4ITStepProcessor::GetProcessNumber()", "Tracking0011",
                        FatalException, "Process Manager is not found.");
            return;
        }

        ActiveOnlyITProcess(pm);
    }
}
// ******************************************************************

void G4ITStepProcessor::ActiveOnlyITProcess(G4ProcessManager* processManager)
{
    G4ProcessVector* processVector = processManager->GetProcessList();

    G4VITProcess* itProcess = 0 ;
    for(int i = 0 ; i < processVector->size() ; i++)
    {
        G4VProcess* base_process = (*processVector)[i];
        itProcess = dynamic_cast<G4VITProcess*>(base_process);

        if(!itProcess)
        {
            int index = processManager->GetProcessIndex(itProcess) ;
            if(index != -1)
                processManager->SetProcessActivation(index, false);
        }
    }
}
// ******************************************************************

void G4ITStepProcessor::SetTrack(G4Track* track)
{
  fpTrack = track ;
  if(fpTrack)
    {
      fpITrack = GetIT(fpTrack) ;
      fpStep = const_cast<G4Step*>(fpTrack -> GetStep());

      if(fpITrack)
        {
          fpTrackingInfo = fpITrack->GetTrackingInfo() ;
        }
      else
        {
            fpTrackingInfo = 0;
            G4cerr << "Track ID : " << fpTrack->GetTrackID() << G4endl;

            G4ExceptionDescription exceptionDescription ("No IT pointer was attached to the track you try to process.");
            G4Exception("G4ITStepProcessor::SetTrack","ITStepProcessor001",
                        FatalErrorInArgument,exceptionDescription);
        }
    }
  else
    {
      fpITrack = 0;
      fpStep = 0 ;
    }
}


void G4ITStepProcessor::GetProcessNumber()
{
#ifdef debug
    G4cout<<"G4ITStepProcessor::GetProcessNumber: is called track="
         <<fTrack<<G4endl;
#endif

    G4ProcessManager* pm= fpTrack->GetDefinition()->GetProcessManager();
    if(!pm)
    {
        G4cerr << "ERROR - G4SteppingManager::GetProcessNumber()" << G4endl
               << "        ProcessManager is NULL for particle = "
               << fpTrack->GetDefinition()->GetParticleName() << ", PDG_code = "
               << fpTrack->GetDefinition()->GetPDGEncoding() << G4endl;
        G4Exception("G4SteppingManager::GetProcessNumber()", "Tracking0011",
                    FatalException, "Process Manager is not found.");
        return;
    }

    // AtRestDoits
    MAXofAtRestLoops =        pm->GetAtRestProcessVector()->entries();
    fpAtRestDoItVector =       pm->GetAtRestProcessVector(typeDoIt);
    fpAtRestGetPhysIntVector = pm->GetAtRestProcessVector(typeGPIL);
#ifdef debug
    G4cout << "G4ITStepProcessor::GetProcessNumber: #ofAtRest="
           << MAXofAtRestLoops << G4endl;
#endif

    // AlongStepDoits
    MAXofAlongStepLoops = pm->GetAlongStepProcessVector()->entries();
    fpAlongStepDoItVector = pm->GetAlongStepProcessVector(typeDoIt);
    fpAlongStepGetPhysIntVector = pm->GetAlongStepProcessVector(typeGPIL);
#ifdef debug
    G4cout << "G4ITStepProcessor::GetProcessNumber:#ofAlongStp="
           << MAXofAlongStepLoops << G4endl;
#endif

    // PostStepDoits
    MAXofPostStepLoops = pm->GetPostStepProcessVector()->entries();
    fpPostStepDoItVector = pm->GetPostStepProcessVector(typeDoIt);
    fpPostStepGetPhysIntVector = pm->GetPostStepProcessVector(typeGPIL);
#ifdef debug
    G4cout << "G4ITStepProcessor::GetProcessNumber: #ofPostStep="
           << MAXofPostStepLoops << G4endl;
#endif

    if (SizeOfSelectedDoItVector<MAXofAtRestLoops    ||
            SizeOfSelectedDoItVector<MAXofAlongStepLoops ||
            SizeOfSelectedDoItVector<MAXofPostStepLoops  )
    {
        G4cerr << "ERROR - G4ITStepProcessor::GetProcessNumber()" << G4endl
               << "        SizeOfSelectedDoItVector= " << SizeOfSelectedDoItVector
               << " ; is smaller then one of MAXofAtRestLoops= "
               << MAXofAtRestLoops << G4endl
               << "        or MAXofAlongStepLoops= " << MAXofAlongStepLoops
               << " or MAXofPostStepLoops= " << MAXofPostStepLoops << G4endl;
        G4Exception("G4ITStepProcessor::GetProcessNumber()",
                    "Tracking0012", FatalException,
                    "The array size is smaller than the actual No of processes.");
    }
}
//______________________________________________________________________________

void G4ITStepProcessor::SetupMembers()
{
    fpSecondary      = fpStep->GetfSecondary();
    fpPreStepPoint   = fpStep->GetPreStepPoint();
    fpPostStepPoint  = fpStep->GetPostStepPoint();

    // Info registered in TrackInformation
    fpPhysicalStep                      = &(fpTrackingInfo->fPhysicalStep) ;
    fpStepStatus                         = &(fpTrackingInfo->fStepStatus) ;
    fpSelectedAtRestDoItVector 		= &(fpTrackingInfo->fSelectedAtRestDoItVector);
    fpSelectedAlongStepDoItVector 	= &(fpTrackingInfo->fSelectedAlongStepDoItVector);
    fpSelectedPostStepDoItVector 	= &(fpTrackingInfo->fSelectedPostStepDoItVector);
    fpTouchableHandle                    = &(fpTrackingInfo->fTouchableHandle);

    GetProcessNumber();
    ResetSecondaries();
}
//______________________________________________________________________________

void G4ITStepProcessor::ResetSecondaries()
{
    // Reset the secondary particles
    fN2ndariesAtRestDoIt    = 0;
    fN2ndariesAlongStepDoIt = 0;
    fN2ndariesPostStepDoIt  = 0;
}
//______________________________________________________________________________

void G4ITStepProcessor::GetAtRestIL()
{
    // Select the rest process which has the shortest time before
    // it is invoked. In rest processes, GPIL()
    // returns the time before a process occurs.
    G4double lifeTime (DBL_MAX), shortestLifeTime (DBL_MAX);

    fAtRestDoItProcTriggered = 0;
    shortestLifeTime = DBL_MAX;

    unsigned int NofInactiveProc=0;

    for( size_t ri=0 ; ri < MAXofAtRestLoops ; ri++ )
    {
        fpCurrentProcess = (G4VITProcess*) (*fpAtRestGetPhysIntVector)[ri];
        if (fpCurrentProcess== 0)
        {
            (*fpSelectedAtRestDoItVector)[ri] = InActivated;
            NofInactiveProc++;
            continue;
        }   // NULL means the process is inactivated by a user on fly.

        fCondition=NotForced;
        lifeTime = fpCurrentProcess->AtRestGPIL( *fpTrack, &fCondition );

        if(fCondition==Forced && fpCurrentProcess)
        {
            (*fpSelectedAtRestDoItVector)[ri] = Forced;
        }
        else
        {
            (*fpSelectedAtRestDoItVector)[ri] = InActivated;
            if(lifeTime < shortestLifeTime )
            {
                shortestLifeTime = lifeTime;
                fAtRestDoItProcTriggered =  G4int(ri);
                (*fpSelectedAtRestDoItVector)[fAtRestDoItProcTriggered] = NotForced;
            }
        }
    }

    fTimeStep = shortestLifeTime ;

    // at least one process is necessary to destroy the particle
    // exit with warning
    if(NofInactiveProc==MAXofAtRestLoops)
    {
        G4cerr << "ERROR - G4ITStepProcessor::InvokeAtRestDoItProcs()" << G4endl
               << "        No AtRestDoIt process is active!" << G4endl;
    }
}
//___________________________________________________________________________

void G4ITStepProcessor::DefinePhysicalStepLength(G4Track* track)
{
    SetTrack(track);
    DoDefinePhysicalStepLength();
}
//______________________________________________________________________________

void G4ITStepProcessor::InitDefineStep()
{

    if(!fpStep)
    {
        SetInitialStep();
    }
    else
    {
        SetupMembers();
        fPreviousStepSize = fpTrack->GetStepLength();
        fPreviousStepTime = fpTrack->GetStep()->GetPostStepPoint()->GetGlobalTime()
                - fpTrack->GetStep()->GetPreStepPoint()->GetGlobalTime() ;

        // Send G4Step information to Hit/Dig if the volume is sensitive
        fpCurrentVolume = fpStep->GetPreStepPoint()->GetPhysicalVolume();
        StepControlFlag =  fpStep->GetControlFlag();
        if( fpCurrentVolume != 0 && StepControlFlag != AvoidHitInvocation)
        {
            fpSensitive = fpStep->GetPreStepPoint()->
                    GetSensitiveDetector();
        }

        // Store last PostStepPoint to PreStepPoint, and swap current and next
        // volume information of G4Track. Reset total energy deposit in one Step.
        fpStep->CopyPostToPreStepPoint();
        fpStep->ResetTotalEnergyDeposit();

        // Switch next touchable in track to current one
        fpTrack->SetTouchableHandle(fpTrack->GetNextTouchableHandle());

        //JA Set the volume before it is used (in DefineStepLength() for User Limit)
        fpCurrentVolume = fpStep->GetPreStepPoint()->GetPhysicalVolume();

        // Reset the step's auxiliary points vector pointer
        fpStep->SetPointerToVectorOfAuxiliaryPoints(0);
    }
}

//______________________________________________________________________________


// ************************************************************************
//	Compute Interaction Length
// ************************************************************************
void G4ITStepProcessor::DoDefinePhysicalStepLength()
{

    InitDefineStep();

    G4TrackStatus trackStatus = fpTrack -> GetTrackStatus()  ;

    if(trackStatus == fStopAndKill)
    {
        return ;
    }

    if(trackStatus == fStopButAlive)
    {
        return GetAtRestIL() ;
    }


    // Find minimum Step length and corresponding time
    // demanded by active disc./cont. processes

    // ReSet the counter etc.
    *fpPhysicalStep  = DBL_MAX;          // Initialize by a huge number
    physIntLength = DBL_MAX;          // Initialize by a huge number

    // GPIL for PostStep
    fPostStepDoItProcTriggered = MAXofPostStepLoops;

    for(size_t np=0; np < MAXofPostStepLoops; np++)
    {
        fpCurrentProcess = (G4VITProcess*) (*fpPostStepGetPhysIntVector)[np];
        if (fpCurrentProcess== 0)
        {
            (*fpSelectedPostStepDoItVector)[np] = InActivated;
            continue;
        }   // NULL means the process is inactivated by a user on fly.

        fCondition=NotForced;
        fpCurrentProcess->SetProcessState(fpTrackingInfo->GetProcessState(fpCurrentProcess->GetProcessID()));

        physIntLength = fpCurrentProcess->
                PostStepGPIL( *fpTrack,
                              fPreviousStepSize,
                              &fCondition );
        fpCurrentProcess->SetProcessState(0);

        switch (fCondition)
        {
        case ExclusivelyForced:
            (*fpSelectedPostStepDoItVector)[np] = ExclusivelyForced;
            *fpStepStatus = fExclusivelyForcedProc;
            fpStep->GetPostStepPoint()
                    ->SetProcessDefinedStep(fpCurrentProcess);
            break;
        case Conditionally:
            (*fpSelectedPostStepDoItVector)[np] = Conditionally;
            break;
        case Forced:
            (*fpSelectedPostStepDoItVector)[np] = Forced;
            break;
        case StronglyForced:
            (*fpSelectedPostStepDoItVector)[np] = StronglyForced;
            break;
        default:
            (*fpSelectedPostStepDoItVector)[np] = InActivated;
            break;
        }

        if (fCondition==ExclusivelyForced)
        {
            for(size_t nrest=np+1; nrest < MAXofPostStepLoops; nrest++)
            {
                (*fpSelectedPostStepDoItVector)[nrest] = InActivated;
            }
            return;  // Please note the 'return' at here !!!
        }
        else
        {
            if(physIntLength < *fpPhysicalStep )
            {
                *fpPhysicalStep = physIntLength;
                *fpStepStatus = fPostStepDoItProc;
                fPostStepDoItProcTriggered = G4int(np);
                fpStep->GetPostStepPoint()
                        ->SetProcessDefinedStep(fpCurrentProcess);
            }
        }
    }

    if (fPostStepDoItProcTriggered<MAXofPostStepLoops)
    {
        if ((*fpSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] ==
                InActivated)
        {
            (*fpSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] =
                    NotForced;
        }
    }

    // GPIL for AlongStep
    proposedSafety = DBL_MAX;
    G4double safetyProposedToAndByProcess = proposedSafety;

    for(size_t kp=0; kp < MAXofAlongStepLoops; kp++)
    {
        fpCurrentProcess = (G4VITProcess*) (*fpAlongStepGetPhysIntVector)[kp];
        if (fpCurrentProcess== 0) continue;
        // NULL means the process is inactivated by a user on fly.

        fpCurrentProcess->SetProcessState(fpTrackingInfo->GetProcessState(fpCurrentProcess->GetProcessID()));
        physIntLength = fpCurrentProcess-> AlongStepGPIL( *fpTrack,
                                                         fPreviousStepSize,
                                                         *fpPhysicalStep,
                                                         safetyProposedToAndByProcess,
                                                         &fGPILSelection );

        if(physIntLength < *fpPhysicalStep)
        {
            *fpPhysicalStep 	= physIntLength;
            // Should save PS and TS in IT

            // Check if the process wants to be the GPIL winner. For example,
            // multi-scattering proposes Step limit, but won't be the winner.
            if(fGPILSelection==CandidateForSelection)
            {
                *fpStepStatus = fAlongStepDoItProc;
                fpStep->GetPostStepPoint()
                        ->SetProcessDefinedStep(fpCurrentProcess);
            }

            // Transportation is assumed to be the last process in the vector
            if(kp == MAXofAlongStepLoops-1)
            {
                fpTransportation = dynamic_cast<G4ITTransportation*>(fpCurrentProcess);

                if(! fpTransportation)
                {
                    G4ExceptionDescription exceptionDescription ;
                    exceptionDescription << "No transportation process found " ;
                    G4Exception("G4ITStepProcessor::DoDefinePhysicalStepLength","G4ITStepProcessor001",
                                FatalErrorInArgument,exceptionDescription);
                }

                fTimeStep 		= fpTransportation->GetInteractionTimeLeft();


                if (fpTrack->GetNextVolume() != 0)
                    *fpStepStatus = fGeomBoundary;
                else
                    *fpStepStatus = fWorldBoundary;
            }
        }
        else
        {
            if(kp == MAXofAlongStepLoops-1)
            {
                fpTransportation = dynamic_cast<G4ITTransportation*>(fpCurrentProcess);

                if(! fpTransportation)
                {
                    G4ExceptionDescription exceptionDescription ;
                    exceptionDescription << "No transportation process found " ;
                    G4Exception("G4ITStepProcessor::DoDefinePhysicalStepLength","G4ITStepProcessor002",
                                FatalErrorInArgument,exceptionDescription);
                }

                fTimeStep 		= fpTransportation->GetInteractionTimeLeft();
            }
        }

        fpCurrentProcess->SetProcessState(0);

        // Make sure to check the safety, even if Step is not limited
        //  by this process.                      J. Apostolakis, June 20, 1998
        //
        if (safetyProposedToAndByProcess < proposedSafety)
            // proposedSafety keeps the smallest value:
            proposedSafety               = safetyProposedToAndByProcess;
        else
            // safetyProposedToAndByProcess always proposes a valid safety:
            safetyProposedToAndByProcess = proposedSafety;

    }
}
//______________________________________________________________________________

void G4ITStepProcessor::Stepping(G4Track* track, const double & timeStep)
{
    CleanProcessor();
    fTimeStep = timeStep ;
    SetTrack(track);
    DoStepping();
}
//______________________________________________________________________________

void G4ITStepProcessor::InitStepping(G4Track* /*track*/)
{
    SetupMembers();

}
//______________________________________________________________________________


// ************************************************************************
//	Stepping
// ************************************************************************
G4StepStatus G4ITStepProcessor::DoStepping()
{

    SetupMembers() ;

    if(!fpAtRestDoItVector  	&&
            !fpAlongStepDoItVector    &&
            !fpPostStepDoItVector)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "No DoIt process found " ;
        G4Exception("G4ITStepProcessor::DoStepping","G4ITStepProcessor003",
                    FatalErrorInArgument,exceptionDescription);
        return fUndefined;
    }
    else if(fpTrack->GetTrackStatus() == fStopAndKill )
    {
        return fUndefined;
    }

    //---------------------------------
    // AtRestStep, AlongStep and PostStep Processes
    //---------------------------------
    else  if( fpTrack->GetTrackStatus() == fStopButAlive )
    {
        if( MAXofAtRestLoops>0 && fpAtRestDoItVector != 0) // second condition to make coverity happy
        {
            //-----------------
            // AtRestStepDoIt
            //-----------------
            InvokeAtRestDoItProcs();
            *fpStepStatus = fAtRestDoItProc;
            fpStep->GetPostStepPoint()->SetStepStatus( *fpStepStatus );

        }
        // Make sure the track is killed
        fpTrack->SetTrackStatus( fStopAndKill );
    }
    else if(fTimeStep > 0.) // Bye, because PostStepIL can return 0 => time =0
    {
        if(fpITrack == 0)
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription
            << " !!! TrackID : "<<  fpTrack->GetTrackID() << G4endl
            << " !!! Track status : "<<  fpTrack->GetTrackStatus() << G4endl
            << " !!! Particle Name : "<< fpTrack -> GetDefinition() -> GetParticleName() << G4endl
            << "No G4ITStepProcessor::fpITrack found" << G4endl;

            G4Exception("G4ITStepProcessor::DoStepping","ITStepProcessor002",
                        FatalErrorInArgument,exceptionDescription);
            return fUndefined; // to make coverity happy
        }

        if(fpITrack->GetTrackingInfo()->IsLeadingStep() == false)
        {
            // In case the track has NOT the minimum step length
            // Given the final step time, the transportation
            // will compute the final position of the particle
            FindTransportationStep();
        }


        // Store the Step length (geometrical length) to G4Step and G4Track
        fpTrack->SetStepLength( *fpPhysicalStep );
        fpStep->SetStepLength( *fpPhysicalStep );

        G4double GeomStepLength = *fpPhysicalStep;

        // Store StepStatus to PostStepPoint
        fpStep->GetPostStepPoint()->SetStepStatus( *fpStepStatus );

        // Invoke AlongStepDoIt
        InvokeAlongStepDoItProcs();

        // Update track by taking into account all changes by AlongStepDoIt
        fpStep->UpdateTrack();

        // Update safety after invocation of all AlongStepDoIts
        endpointSafOrigin= fpPostStepPoint->GetPosition();

        endpointSafety=  std::max( proposedSafety - GeomStepLength, kCarTolerance);

        fpStep->GetPostStepPoint()->SetSafety( endpointSafety );

        if(GetIT(fpTrack)->GetTrackingInfo()->IsLeadingStep())
        {
            // Invoke PostStepDoIt including G4ITTransportation::PSDI
            InvokePostStepDoItProcs();
        }
        else
        {
            // Only invoke transportation
            InvokeTransportationProc();
        }
    }

    //-------
    // Finale
    //-------

    // Update 'TrackLength' and remeber the Step length of the current Step
    fpTrack->AddTrackLength(fpStep->GetStepLength());
    fpTrack->IncrementCurrentStepNumber();

    // Send G4Step information to Hit/Dig if the volume is sensitive
    fpCurrentVolume = fpStep->GetPreStepPoint()->GetPhysicalVolume();
    StepControlFlag =  fpStep->GetControlFlag();
/***
    if( fpCurrentVolume != 0 && StepControlFlag != AvoidHitInvocation)
    {
        fpSensitive = fpStep->GetPreStepPoint()->
                GetSensitiveDetector();
        if( fpSensitive != 0 )
        {
            fpSensitive->Hit(fpStep);
        }
    }

     User intervention process.
    if( fpUserSteppingAction != 0 )
    {
        fpUserSteppingAction->UserSteppingAction(fpStep);
    }
    G4UserSteppingAction* regionalAction
            = fpStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetRegion()
            ->GetRegionalSteppingAction();
    if( regionalAction ) regionalAction->UserSteppingAction(fpStep);
***/
    fpTrackingManager->AppendTrajectory(fpTrack,fpStep);
    // Stepping process finish. Return the value of the StepStatus.
    return *fpStepStatus;
}

///////////////////////////////////////////////////////////
void G4ITStepProcessor::SetInitialStep()
///////////////////////////////////////////////////////////
{
    // DEBUG
//    G4cout << "SetInitialStep for : " << fpITrack-> GetName() << G4endl;
    fpStep = new G4Step();
    fpTrack->SetStep(fpStep);
    fpSecondary = fpStep->NewSecondaryVector();

    SetupMembers();

    //________________________________________________________
    // Initialize geometry
    if ( ! *fpTouchableHandle )
    {
        G4ThreeVector direction= fpTrack->GetMomentumDirection();
        fpNavigator->LocateGlobalPointAndSetup( fpTrack->GetPosition(),
                                               &direction, false, false );
        *fpTouchableHandle = fpNavigator->CreateTouchableHistory();

        fpTrack->SetTouchableHandle( *fpTouchableHandle );
        fpTrack->SetNextTouchableHandle( *fpTouchableHandle );
    }
    else
    {
        fpTrack->SetNextTouchableHandle( *fpTouchableHandle );
        G4VPhysicalVolume* oldTopVolume= fpTrack->GetTouchableHandle()->GetVolume();
        G4VPhysicalVolume* newTopVolume=
                fpNavigator->ResetHierarchyAndLocate( fpTrack->GetPosition(),
                                                     fpTrack->GetMomentumDirection(),
                                                     *((G4TouchableHistory*)fpTrack->GetTouchableHandle()()) );
        if(newTopVolume != oldTopVolume || oldTopVolume->GetRegularStructureId() == 1 )
        {
            *fpTouchableHandle = fpNavigator->CreateTouchableHistory();
            fpTrack->SetTouchableHandle( *fpTouchableHandle );
            fpTrack->SetNextTouchableHandle( *fpTouchableHandle );
        }
    }

    fpCurrentVolume = (*fpTouchableHandle)->GetVolume();

    //________________________________________________________
    // If the primary track has 'Suspend' or 'PostponeToNextEvent' state,
    // set the track state to 'Alive'.
    if( (fpTrack->GetTrackStatus()==fSuspend) ||
            (fpTrack->GetTrackStatus()==fPostponeToNextEvent) )
    {
        fpTrack->SetTrackStatus(fAlive);
    }

    // If the primary track has 'zero' kinetic energy, set the track
    // state to 'StopButAlive'.
    if(fpTrack->GetKineticEnergy() <= 0.0)
    {
        fpTrack->SetTrackStatus( fStopButAlive );
    }
    //________________________________________________________
    // Set vertex information of G4Track at here
    if ( fpTrack->GetCurrentStepNumber() == 0 )
    {
        fpTrack->SetVertexPosition( fpTrack->GetPosition() );
        fpTrack->SetVertexMomentumDirection( fpTrack->GetMomentumDirection() );
        fpTrack->SetVertexKineticEnergy( fpTrack->GetKineticEnergy() );
        fpTrack->SetLogicalVolumeAtVertex( fpTrack->GetVolume()->GetLogicalVolume() );
    }
    //________________________________________________________
    // If track is already outside the world boundary, kill it
    if( fpCurrentVolume==0 )
    {
        // If the track is a primary, stop processing
        if(fpTrack->GetParentID()==0)
        {
            G4cerr << "ERROR - G4ITStepProcessor::SetInitialStep()" << G4endl
                   << "        Primary particle starting at - "
                   << fpTrack->GetPosition()
                   << " - is outside of the world volume." << G4endl;
            G4Exception("G4ITStepProcessor::SetInitialStep()", "ITStepProcessor003",
                        FatalException, "Primary vertex outside of the world!");
        }

        fpTrack->SetTrackStatus( fStopAndKill );
        G4cout << "WARNING - G4ITStepProcessor::SetInitialStep()" << G4endl
               << "          Initial track position is outside world! - "
               << fpTrack->GetPosition() << G4endl;
    }
    else{
        // Initial set up for attribues of 'Step'
        fpStep->InitializeStep( fpTrack );
    }


    if( fpTrack->GetTrackStatus() == fStopAndKill ) return ;

    fpTrackingManager->StartTracking(fpTrack);

    *fpStepStatus = fUndefined;

}
//______________________________________________________________________________

void G4ITStepProcessor::FindTransportationStep()
{
    double physicalStep = DBL_MAX ;

    fpTransportation = dynamic_cast<G4ITTransportation*>
            ((*fpAlongStepGetPhysIntVector)[MAXofAlongStepLoops-1]);

    if(!fpTrack)
    {        
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription
        << "No G4ITStepProcessor::fpTrack found";
        G4Exception("G4ITStepProcessor::FindTransportationStep","ITStepProcessor004",
                    FatalErrorInArgument,exceptionDescription);
        return;

    }
    if(!fpITrack)
    {        
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription
        << "No G4ITStepProcessor::fITrack" ;
        G4Exception("G4ITStepProcessor::FindTransportationStep","ITStepProcessor005",
                    FatalErrorInArgument,exceptionDescription);
        return;
    }
    if(!(fpITrack->GetTrack()))
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription
        << "No G4ITStepProcessor::fITrack->GetTrack()" ;
        G4Exception("G4ITStepProcessor::FindTransportationStep","ITStepProcessor006",
                    FatalErrorInArgument,exceptionDescription);
        return;
    }

    if(fpTransportation)
    {
        fpTransportation->SetProcessState(fpTrackingInfo->GetProcessState(fpTransportation->GetProcessID()));
        fpTransportation -> ComputeStep(*fpTrack, *fpStep, fTimeStep, physicalStep) ;
        fpTransportation->SetProcessState(0);
    }

    if(physicalStep >= DBL_MAX)
    {
        fpTrack -> SetTrackStatus(fStopAndKill) ;
        return ;
    }

    *fpPhysicalStep = physicalStep ;
}
//______________________________________________________________________________

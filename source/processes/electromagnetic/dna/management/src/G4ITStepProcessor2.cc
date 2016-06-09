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
#include "G4LossTableManager.hh"
#include "G4EnergyLossTables.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VITProcess.hh"
#include "G4TrackingInformation.hh"


void G4ITStepProcessor::DealWithSecondaries(G4int& counter)
{
    // Now Store the secondaries from ParticleChange to SecondaryList
    G4Track* tempSecondaryTrack;

    for(G4int DSecLoop=0 ;
        DSecLoop<fpParticleChange->GetNumberOfSecondaries() ;
        DSecLoop++)
    {
        tempSecondaryTrack = fpParticleChange->GetSecondary(DSecLoop);

        if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
        {
            ApplyProductionCut(tempSecondaryTrack);
        }

        // Set parentID
        tempSecondaryTrack->SetParentID( fpTrack->GetTrackID() );

        // Set the process pointer which created this track
        tempSecondaryTrack->SetCreatorProcess( fpCurrentProcess );

        // If this 2ndry particle has 'zero' kinetic energy, make sure
        // it invokes a rest process at the beginning of the tracking
        if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN)
        {
            G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()->GetProcessManager();
            if (pm->GetAtRestProcessVector()->entries()>0){
              tempSecondaryTrack->SetTrackStatus( fStopButAlive );
              fpSecondary->push_back( tempSecondaryTrack );
              fN2ndariesAtRestDoIt++;
            } else {
              delete tempSecondaryTrack;
            }
        }
        else
        {
            fpSecondary->push_back( tempSecondaryTrack );
            counter++;
        }
    } //end of loop on secondary
}

//////////////////////////////////////////////////////
void G4ITStepProcessor::InvokeAtRestDoItProcs()
//////////////////////////////////////////////////////
{
    fpStep->SetStepLength( 0. );  //the particle has stopped
    fpTrack->SetStepLength( 0. );

// invoke selected process
    for(size_t np=0; np < MAXofAtRestLoops; np++)
    {
        //
        // Note: DoItVector has inverse order against GetPhysIntVector
        //       and SelectedAtRestDoItVector.
        //
        if( (*fpSelectedAtRestDoItVector)[MAXofAtRestLoops-np-1] != InActivated)
        {
            fpCurrentProcess = (G4VITProcess*) (*fpAtRestDoItVector)[np];

            fpCurrentProcess->SetProcessState(
                        fpTrackingInfo->GetProcessState(fpCurrentProcess->GetProcessID()));
            fpParticleChange
            = fpCurrentProcess->AtRestDoIt( *fpTrack, *fpStep);
            fpCurrentProcess->SetProcessState(0);

            // Set the current process as a process which defined this Step length
            fpStep->GetPostStepPoint()
            ->SetProcessDefinedStep(fpCurrentProcess);

            // Update Step
            fpParticleChange->UpdateStepForAtRest(fpStep);

            // Now Store the secondaries from ParticleChange to SecondaryList
            DealWithSecondaries(fN2ndariesAtRestDoIt) ;

            // clear ParticleChange
            fpParticleChange->Clear();

        } //if(fSelectedAtRestDoItVector[np] != InActivated){
    } //for(size_t np=0; np < MAXofAtRestLoops; np++){
    fpStep->UpdateTrack();

    fpTrack->SetTrackStatus( fStopAndKill );
}
//______________________________________________________________________________

void G4ITStepProcessor::InvokeAlongStepDoItProcs()
{

// If the current Step is defined by a 'ExclusivelyForced'
// PostStepDoIt, then don't invoke any AlongStepDoIt
    if(*fpStepStatus == fExclusivelyForcedProc)
    {
        return;               // Take note 'return' at here !!!
    }

// Invoke the all active continuous processes
    for( size_t ci=0 ; ci<MAXofAlongStepLoops ; ci++ )
    {
        fpCurrentProcess = (G4VITProcess*) (*fpAlongStepDoItVector)[ci];
        if (fpCurrentProcess== 0) continue;
        // NULL means the process is inactivated by a user on fly.

        fpCurrentProcess->SetProcessState(fpTrackingInfo->GetProcessState(fpCurrentProcess->GetProcessID()));
        fpParticleChange
        = fpCurrentProcess->AlongStepDoIt( *fpTrack, *fpStep );
        fpCurrentProcess->SetProcessState(0);
        // Update the PostStepPoint of Step according to ParticleChange
        fpParticleChange->UpdateStepForAlongStep(fpStep);

        // Now Store the secondaries from ParticleChange to SecondaryList
        DealWithSecondaries(fN2ndariesAlongStepDoIt) ;

        // Set the track status according to what the process defined
        // if kinetic energy >0, otherwise set  fStopButAlive
        fpTrack->SetTrackStatus( fpParticleChange->GetTrackStatus() );

        // clear ParticleChange
        fpParticleChange->Clear();
    }

    fpStep->UpdateTrack();
    G4TrackStatus fNewStatus = fpTrack->GetTrackStatus();

    if ( fNewStatus == fAlive && fpTrack->GetKineticEnergy() <= DBL_MIN )
    {
        G4cout << "G4ITStepProcessor::InvokeAlongStepDoItProcs : Track will be killed" << G4endl;
        if(MAXofAtRestLoops>0) fNewStatus = fStopButAlive;
        else                   fNewStatus = fStopAndKill;
        fpTrack->SetTrackStatus( fNewStatus );
    }

}

////////////////////////////////////////////////////////
void G4ITStepProcessor::InvokePostStepDoItProcs()
////////////////////////////////////////////////////////
{

// Invoke the specified discrete processes
    for(size_t np=0; np < MAXofPostStepLoops; np++)
    {
        //
        // Note: DoItVector has inverse order against GetPhysIntVector
        //       and SelectedPostStepDoItVector.
        //
        G4int Cond = (*fpSelectedPostStepDoItVector)[MAXofPostStepLoops-np-1];
        if(Cond != InActivated)
        {
            if( ((Cond == NotForced) && (*fpStepStatus == fPostStepDoItProc)) ||
                    ((Cond == Forced) && (*fpStepStatus != fExclusivelyForcedProc)) ||
                    ((Cond == Conditionally) && (*fpStepStatus == fAlongStepDoItProc)) ||
                    ((Cond == ExclusivelyForced) && (*fpStepStatus == fExclusivelyForcedProc)) ||
                    ((Cond == StronglyForced) )
              )
            {

                InvokePSDIP(np);
            }
        } //if(*fSelectedPostStepDoItVector(np)........

        // Exit from PostStepLoop if the track has been killed,
        // but extra treatment for processes with Strongly Forced flag
        if(fpTrack->GetTrackStatus() == fStopAndKill)
        {
            for(size_t np1=np+1; np1 < MAXofPostStepLoops; np1++)
            {
                G4int Cond2 = (*fpSelectedPostStepDoItVector)[MAXofPostStepLoops-np1-1];
                if (Cond2 == StronglyForced)
                {
                    InvokePSDIP(np1);
                }
            }
            break;
        }
    } //for(size_t np=0; np < MAXofPostStepLoops; np++){
}

////////////////////////////////////////////////////////
void G4ITStepProcessor::InvokePSDIP(size_t np)
////////////////////////////////////////////////////////
{
    fpCurrentProcess = (G4VITProcess*) (*fpPostStepDoItVector)[np];

    fpCurrentProcess->SetProcessState(fpTrackingInfo->GetProcessState(fpCurrentProcess->GetProcessID()));
    fpParticleChange
    = fpCurrentProcess->PostStepDoIt( *fpTrack, *fpStep);
    fpCurrentProcess->SetProcessState(0);

    // Update PostStepPoint of Step according to ParticleChange
    fpParticleChange->UpdateStepForPostStep(fpStep);

    // Update G4Track according to ParticleChange after each PostStepDoIt
    fpStep->UpdateTrack();

    // Update safety after each invocation of PostStepDoIts
    fpStep->GetPostStepPoint()->SetSafety( CalculateSafety() );

    // Now Store the secondaries from ParticleChange to SecondaryList
    DealWithSecondaries(fN2ndariesPostStepDoIt) ;

    // Set the track status according to what the process defined
    fpTrack->SetTrackStatus( fpParticleChange->GetTrackStatus() );

    // clear ParticleChange
    fpParticleChange->Clear();
}

////////////////////////////////////////////////////////
void G4ITStepProcessor::InvokeTransportationProc()
////////////////////////////////////////////////////////
{
    // Invoke the specified discrete processes
    for(size_t np=0; np < MAXofPostStepLoops; np++)
    {
        //
        // Note: DoItVector has inverse order against GetPhysIntVector
        //       and SelectedPostStepDoItVector.
        //
        G4int Cond = (*fpSelectedPostStepDoItVector)[MAXofPostStepLoops-np-1];
        if(Cond != InActivated)
        {
            if(
                ((Cond == Forced) && (*fpStepStatus != fExclusivelyForcedProc)) ||
                ((Cond == Conditionally) && (*fpStepStatus == fAlongStepDoItProc)) ||
                ((Cond == ExclusivelyForced) && (*fpStepStatus == fExclusivelyForcedProc)) ||
                ((Cond == StronglyForced) )
              )
            {

                InvokePSDIP(np);
            }
        } //if(*fSelectedPostStepDoItVector(np)........

        // Exit from PostStepLoop if the track has been killed,
        // but extra treatment for processes with Strongly Forced flag
        if(fpTrack->GetTrackStatus() == fStopAndKill)
        {
            for(size_t np1=np+1; np1 < MAXofPostStepLoops; np1++)
            {
                G4int Cond2 = (*fpSelectedPostStepDoItVector)[MAXofPostStepLoops-np1-1];
                if (Cond2 == StronglyForced)
                {
                    InvokePSDIP(np1);
                }
            }
            break;
        }
    }
}

///////////////////////////////////////////////////////////////
void G4ITStepProcessor::ApplyProductionCut(G4Track* aSecondary)
///////////////////////////////////////////////////////////////
{
    G4bool tBelowCutEnergyAndSafety = false;
    G4int tPtclIdx
    = G4ProductionCuts::GetIndex(aSecondary->GetDefinition());
    if (tPtclIdx<0)
    {
        return;
    }
    G4ProductionCutsTable* tCutsTbl
    = G4ProductionCutsTable::GetProductionCutsTable();
    G4int tCoupleIdx
    = tCutsTbl->GetCoupleIndex(fpPreStepPoint->GetMaterialCutsCouple());
    G4double tProdThreshold
    = (*(tCutsTbl->GetEnergyCutsVector(tPtclIdx)))[tCoupleIdx];
    if( aSecondary->GetKineticEnergy()<tProdThreshold )
    {
        tBelowCutEnergyAndSafety = true;
        if(std::abs(aSecondary->GetDynamicParticle()->GetCharge()) > DBL_MIN)
        {
            G4double currentRange
            = G4LossTableManager::Instance()->GetRange(aSecondary->GetDefinition(),
                    aSecondary->GetKineticEnergy(),
                    fpPreStepPoint->GetMaterialCutsCouple());
            tBelowCutEnergyAndSafety = (currentRange < CalculateSafety() );
        }
    }

    if( tBelowCutEnergyAndSafety )
    {
        if( !(aSecondary->IsGoodForTracking()) )
        {
            // Add kinetic energy to the total energy deposit
            fpStep->AddTotalEnergyDeposit(
                aSecondary->GetKineticEnergy() );
            aSecondary->SetKineticEnergy(0.0);
        }
    }
}

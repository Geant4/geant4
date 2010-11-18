/*
 * =============================================================================
 *
 *       Filename:  CexmcStudiedProcess.cc
 *
 *    Description:  studied process in the target
 *
 *        Version:  1.0
 *        Created:  14.09.2010 18:37:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <G4VParticleChange.hh>
#include "CexmcStudiedProcess.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcIncidentParticleTrackInfo.hh"
#include "CexmcCommon.hh"


CexmcStudiedProcess::CexmcStudiedProcess( CexmcPhysicsManager *  physicsManager,
                                          G4ProcessType  processType ) :
    G4WrapperProcess( CexmcStudiedProcessFirstName, processType ),
    physicsManager( physicsManager )
{
}


G4double  CexmcStudiedProcess::PostStepGetPhysicalInteractionLength(
            const G4Track &  track, G4double , G4ForceCondition *  condition )
{
    *condition = NotForced;

    if ( ! physicsManager->IsStudiedProcessAllowed() )
        return CexmcDblMax;

    CexmcTrackInfo *  trackInfo( static_cast< CexmcTrackInfo * >(
                                                track.GetUserInformation() ) );

    if ( ! trackInfo ||
         trackInfo->GetTypeInfo() != CexmcIncidentParticleTrackType )
        return CexmcDblMax;

    CexmcIncidentParticleTrackInfo *  theTrackInfo(
                static_cast< CexmcIncidentParticleTrackInfo * >( trackInfo ) );

    if ( ! theTrackInfo->IsStudiedProcessActivated() )
        return CexmcDblMax;

    return theTrackInfo->GetFinalTrackLengthInTarget() -
            theTrackInfo->GetCurrentTrackLengthInTarget();
}


G4VParticleChange *  CexmcStudiedProcess::PostStepDoIt( const G4Track &  track,
                                                        const G4Step &  step )
{
    G4VParticleChange *  particleChange( pRegProcess->PostStepDoIt( track,
                                                                    step ) );

    if ( particleChange && particleChange->GetTrackStatus() == fStopAndKill )
        physicsManager->IncrementNumberOfTriggeredStudiedInteractions();

    return particleChange;
}


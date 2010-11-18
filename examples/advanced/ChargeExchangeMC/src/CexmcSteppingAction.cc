/*
 * ============================================================================
 *
 *       Filename:  CexmcSteppingAction.cc
 *
 *    Description:  stepping action
 *
 *        Version:  1.0
 *        Created:  27.10.2009 16:03:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4Step.hh>
#include <G4Track.hh>
#include <G4ParticleDefinition.hh>
#include <G4VTouchable.hh>
#include <G4NavigationHistory.hh>
#include <G4AffineTransform.hh>
#include <G4UnitsTable.hh>
#include <G4RunManager.hh>
#include "CexmcSteppingAction.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcSetup.hh"
#include "CexmcIncidentParticleTrackInfo.hh"
#include "CexmcCommon.hh"


CexmcSteppingAction::CexmcSteppingAction(
                                    CexmcPhysicsManager *  physicsManager ) :
    physicsManager( physicsManager ), targetVolume( NULL )
{
    G4RunManager *      runManager( G4RunManager::GetRunManager() );
    const CexmcSetup *  setup( static_cast< const CexmcSetup * >(
                                runManager->GetUserDetectorConstruction() ) );
    targetVolume = setup->GetVolume( CexmcSetup::Target );
}


void  CexmcSteppingAction::UserSteppingAction( const G4Step *  step )
{
    G4Track *         track( step->GetTrack() );
    CexmcTrackInfo *  trackInfo( static_cast< CexmcTrackInfo * >(
                                                track->GetUserInformation() ) );

    if ( ! trackInfo ||
         trackInfo->GetTypeInfo() != CexmcIncidentParticleTrackType )
        return;

    CexmcIncidentParticleTrackInfo *  theTrackInfo(
                static_cast< CexmcIncidentParticleTrackInfo * >( trackInfo ) );

    G4StepPoint *         postStepPoint( step->GetPostStepPoint() );
    G4StepStatus          stepStatus( postStepPoint->GetStepStatus() );
    const G4VTouchable *  touchable( postStepPoint->GetTouchable() );
    G4VPhysicalVolume *   volume( touchable->GetVolume() );

    if ( volume && volume->GetLogicalVolume() == targetVolume )
    {
        if ( ! theTrackInfo->IsStudiedProcessActivated() )
        {
            physicsManager->ResampleTrackLengthInTarget( track, postStepPoint );
            theTrackInfo->ActivateStudiedProcess();
        }

        if ( stepStatus != fGeomBoundary )
        {
            if ( theTrackInfo->NeedsTrackLengthResampling() )
                physicsManager->ResampleTrackLengthInTarget(
                                                        track, postStepPoint );
            else
                theTrackInfo->AddTrackLengthInTarget( step->GetStepLength() );
        }
    }

    G4StepPoint *  preStepPoint( step->GetPreStepPoint() );
    touchable = preStepPoint->GetTouchable();
    volume = touchable->GetVolume();

    if ( volume && volume->GetLogicalVolume() == targetVolume )
    {
        if ( stepStatus == fGeomBoundary )
            theTrackInfo->ActivateStudiedProcess( false );
    }
}


/*
 * ============================================================================
 *
 *       Filename:  CexmcTrackingAction.cc
 *
 *    Description:  tracking action
 *
 *        Version:  1.0
 *        Created:  22.11.2009 18:22:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4ParticleDefinition.hh>
#include <G4VProcess.hh>
#include <G4Track.hh>
#include <G4RunManager.hh>
#include "CexmcTrackingAction.hh"
#include "CexmcTrackInfo.hh"
#include "CexmcIncidentParticleTrackInfo.hh"
#include "CexmcProductionModel.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcSetup.hh"
#include "CexmcException.hh"
#include "CexmcCommon.hh"


CexmcTrackingAction::CexmcTrackingAction(
                                    CexmcPhysicsManager *  physicsManager ) :
    physicsManager( physicsManager ), targetVolume( NULL ),
    outputParticleTrackId( CexmcInvalidTrackId ),
    outputParticleDecayProductCopyNumber( 0 ), incidentParticle( NULL ),
    outputParticle( NULL ), nucleusOutputParticle( NULL )
{
    CexmcProductionModel *  productionModel(
                                    physicsManager->GetProductionModel() );
    if ( ! productionModel )
        throw CexmcException( CexmcWeirdException );

    incidentParticle = productionModel->GetIncidentParticle();
    outputParticle = productionModel->GetOutputParticle();
    nucleusOutputParticle = productionModel->GetNucleusOutputParticle();

    if ( ! incidentParticle || ! outputParticle || ! nucleusOutputParticle )
        throw CexmcException( CexmcIncompleteProductionModel );

    G4RunManager *      runManager( G4RunManager::GetRunManager() );
    const CexmcSetup *  setup( static_cast< const CexmcSetup * >(
                                runManager->GetUserDetectorConstruction() ) );
    targetVolume = setup->GetVolume( CexmcSetup::Target );
}


void  CexmcTrackingAction::PreUserTrackingAction( const G4Track *  track )
{
    CexmcTrackInfo *  trackInfo( static_cast< CexmcTrackInfo * >(
                                                track->GetUserInformation() ) );

    if ( trackInfo )
        return;

    G4Track *  theTrack( const_cast< G4Track * >( track ) );

    do
    {
        if ( track->GetParentID() == 0 )
        {
            if ( *track->GetDefinition() == *incidentParticle )
            {
                trackInfo = new CexmcIncidentParticleTrackInfo(
                                                    CexmcBeamParticleTrack );
                theTrack->SetUserInformation( trackInfo );
                SetupIncidentParticleTrackInfo( track );
            }
            else
            {
                trackInfo = new CexmcTrackInfo( CexmcBeamParticleTrack );
            }
            break;
        }

        if ( track->GetCreatorProcess()->GetProcessName() ==
             CexmcStudiedProcessFullName )
        {
            do
            {
                if ( *track->GetDefinition() == *outputParticle )
                {
                    outputParticleTrackId = track->GetTrackID();
                    trackInfo = new CexmcTrackInfo( CexmcOutputParticleTrack );
                    break;
                }
                if ( *track->GetDefinition() == *nucleusOutputParticle )
                {
                    trackInfo = new CexmcTrackInfo( CexmcNucleusParticleTrack );
                    break;
                }
            } while ( false );
            break;
        }

        if ( track->GetParentID() == outputParticleTrackId )
        {
            trackInfo = new CexmcTrackInfo(
                                    CexmcOutputParticleDecayProductTrack,
                                    outputParticleDecayProductCopyNumber++ );
            break;
        }

        if ( *track->GetDefinition() == *incidentParticle )
        {
            if ( physicsManager->OnlyBeamParticleCanTriggerStudiedProcess() )
                break;
            trackInfo = new CexmcIncidentParticleTrackInfo( CexmcInsipidTrack );
            theTrack->SetUserInformation( trackInfo );
            SetupIncidentParticleTrackInfo( track );
            break;
        }
    } while ( false );

    if ( ! trackInfo )
        return;

    if ( ! track->GetUserInformation() )
        theTrack->SetUserInformation( trackInfo );
}


void  CexmcTrackingAction::SetupIncidentParticleTrackInfo(
                                                    const G4Track *  track )
{
    CexmcIncidentParticleTrackInfo *  trackInfo(
                    static_cast< CexmcIncidentParticleTrackInfo * >(
                                                track->GetUserInformation() ) );

    if ( ! trackInfo )
        return;

    G4VPhysicalVolume *  volume( track->GetVolume() );

    if ( volume && volume->GetLogicalVolume() == targetVolume )
    {
        physicsManager->ResampleTrackLengthInTarget( track );
        trackInfo->ActivateStudiedProcess();
    }
}


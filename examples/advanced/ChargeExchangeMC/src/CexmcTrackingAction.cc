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
                                    CexmcPhysicsManager *  physicsManager_ ) :
    physicsManager( physicsManager_ ), targetVolume( NULL ),
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


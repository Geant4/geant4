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
                                    CexmcPhysicsManager *  physicsManager_ ) :
    physicsManager( physicsManager_ ), targetVolume( NULL )
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


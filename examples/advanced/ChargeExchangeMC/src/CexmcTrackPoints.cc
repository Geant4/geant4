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
 *       Filename:  CexmcTrackPoints.cc
 *
 *    Description:  track points collection
 *
 *        Version:  1.0
 *        Created:  16.11.2009 13:03:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4HCofThisEvent.hh>
#include <G4AffineTransform.hh>
#include <G4UnitsTable.hh>
#include "CexmcTrackPoints.hh"
#include "CexmcTrackInfo.hh"


CexmcTrackPoints::CexmcTrackPoints( const G4String &  name ) :
    CexmcPrimitiveScorer( name ), hcId( -1 )
{
}


G4int  CexmcTrackPoints::GetTrackId( G4Step *  step )
{
    /* NB: trackInfo for these tracks shall always exist as far as
     * track points filter decline tracks without it */
    CexmcTrackInfo *  trackInfo( static_cast< CexmcTrackInfo * >(
                                    step->GetTrack()->GetUserInformation() ) );

    G4int  ret( trackInfo->GetTrackType() );

    if ( ret == CexmcOutputParticleDecayProductTrack )
        ret += trackInfo->GetCopyNumber();

    return ret;
}


G4int  CexmcTrackPoints::GetIndex( G4Step *  step )
{
    return GetTrackId( step );
}


G4bool  CexmcTrackPoints::ProcessHits( G4Step *  step, G4TouchableHistory * )
{
    G4int  index( GetIndex( step ) );

    if ( ( *eventMap )[ index ] )
        return false;

    G4Track *               track( step->GetTrack() );
    G4ParticleDefinition *  particle( track->GetDefinition() );
    CexmcTrackType          trackType( CexmcInsipidTrack );

    G4StepPoint *           preStepPoint( step->GetPreStepPoint() );
    G4ThreeVector           position( preStepPoint->GetPosition() );
    G4ThreeVector           direction( preStepPoint->GetMomentumDirection() );

    const G4AffineTransform &  transform( preStepPoint->GetTouchable()->
                                          GetHistory()->GetTopTransform() );

    CexmcTrackInfo *  trackInfo( static_cast< CexmcTrackInfo * >(
                                                track->GetUserInformation() ) );
    trackType = trackInfo->GetTrackType();

    CexmcTrackPointInfo  trackPointInfo( transform.TransformPoint( position ),
                                         position,
                                         transform.TransformAxis( direction ),
                                         direction,
                                         preStepPoint->GetMomentum().mag(),
                                         particle, track->GetTrackID(),
                                         trackType );

    eventMap->set( index, trackPointInfo );

    return true; 
}


void  CexmcTrackPoints::Initialize( G4HCofThisEvent *  hcOfEvent )
{
    eventMap = new CexmcTrackPointsCollection( detector->GetName(),
                                               primitiveName );

    if ( hcId < 0 )
        hcId = GetCollectionID( 0 );

    hcOfEvent->AddHitsCollection( hcId, eventMap );
}


void  CexmcTrackPoints::EndOfEvent( G4HCofThisEvent * )
{
    if ( GetVerboseLevel() > 0 )
        PrintAll();
}


void  CexmcTrackPoints::clear( void )
{
    eventMap->clear();
}


void  CexmcTrackPoints::DrawAll( void )
{
}


void  CexmcTrackPoints::PrintAll( void )
{
    G4int   nmbOfEntries( eventMap->entries() );

    if ( nmbOfEntries == 0 )
        return;

    PrintHeader( nmbOfEntries );

    for ( CexmcTrackPointsCollectionData::iterator
                         itr( eventMap->GetMap()->begin() );
                                     itr != eventMap->GetMap()->end(); ++itr )
    {
        G4cout << "       track id " << itr->first << G4endl;
        G4cout << "         , position: " <<
                G4BestUnit( itr->second->positionLocal, "Length" ) << G4endl;
        G4cout << "         , direction: " <<
                itr->second->directionLocal << G4endl;
        G4cout << "         , momentum: " <<
                G4BestUnit( itr->second->momentumAmp, "Energy" ) << G4endl;
        G4cout << "         , particle: " <<
                itr->second->particle->GetParticleName() << G4endl;
    }
}


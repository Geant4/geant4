/*
 * ============================================================================
 *
 *       Filename:  CexmcTrackPointsInLeftRightSet.cc
 *
 *    Description:  track points in left/right detector sets
 *                  (e.g. veto counters and calorimeters)
 *
 *        Version:  1.0
 *        Created:  22.11.2009 21:15:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4VTouchable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4UnitsTable.hh>
#include "CexmcTrackPointsInLeftRightSet.hh"
#include "CexmcSetup.hh"


G4int  CexmcTrackPointsInLeftRightSet::leftRightBitsOffset( 24 );


CexmcTrackPointsInLeftRightSet::CexmcTrackPointsInLeftRightSet(
                        const G4String &  name, const CexmcSetup *  setup ) :
    CexmcTrackPoints( name ), setup( setup )
{
}


G4int  CexmcTrackPointsInLeftRightSet::GetIndex( G4Step *  step )
{
    G4int                ret( GetTrackId( step ) );
    G4StepPoint *        preStep( step->GetPreStepPoint() );
    G4VPhysicalVolume *  pVolume( preStep->GetPhysicalVolume() );

    if ( setup->IsRightDetector( pVolume ) )
        ret |= 1 << leftRightBitsOffset;

    return ret;
}


void  CexmcTrackPointsInLeftRightSet::PrintAll( void )
{
    G4int   nmbOfEntries( eventMap->entries() );

    if ( nmbOfEntries == 0 )
        return;

    G4cout << " --- MultiFunctionalDet " << detector->GetName() << G4endl;
    G4cout << "     PrimitiveScorer " << GetName() << G4endl;
    G4cout << "     Number of entries " << nmbOfEntries << G4endl;

    for( std::map< G4int, CexmcTrackPointInfo* >::iterator
                                     itr( eventMap->GetMap()->begin() );
         itr != eventMap->GetMap()->end(); ++itr )
    {
        G4bool  isRightDetector( itr->first >> leftRightBitsOffset );
        const G4String  detectorSide( isRightDetector ? "right" : "left" );
        G4int   trackId( itr->first &
                            ( ( 1 << ( leftRightBitsOffset - 1 ) ) |
                              ( ( 1 << ( leftRightBitsOffset - 1 ) ) - 1 ) ) );
        G4cout << "       " << detectorSide << " detector" << G4endl;
        G4cout << "         , track id " << trackId << G4endl;
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


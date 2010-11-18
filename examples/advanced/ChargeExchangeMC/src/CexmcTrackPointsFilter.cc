/*
 * ============================================================================
 *
 *       Filename:  CexmcTrackPointsFilter.cc
 *
 *    Description:  track points of interest
 *
 *        Version:  1.0
 *        Created:  16.11.2009 22:29:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4String.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>
#include "CexmcTrackPointsFilter.hh"
#include "CexmcTrackInfo.hh"
#include "CexmcCommon.hh"


CexmcTrackPointsFilter::CexmcTrackPointsFilter( const G4String &  name ) :
    G4VSDFilter( name )
{
}


G4bool  CexmcTrackPointsFilter::Accept( const G4Step *  step ) const
{
    G4Track *  track( step->GetTrack() );
    CexmcTrackInfo *  trackInfo( static_cast< CexmcTrackInfo * >(
                                                track->GetUserInformation() ) );

    if ( ! trackInfo )
        return false;

    return trackInfo->GetTrackType() != CexmcInsipidTrack;
}


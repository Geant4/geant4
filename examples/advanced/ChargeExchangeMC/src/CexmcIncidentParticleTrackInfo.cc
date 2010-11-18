/*
 * =============================================================================
 *
 *       Filename:  CexmcIncidentParticleTrackInfo.cc
 *
 *    Description:  incident particle track info
 *
 *        Version:  1.0
 *        Created:  18.05.2010 13:18:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include "CexmcIncidentParticleTrackInfo.hh"


CexmcIncidentParticleTrackInfo::CexmcIncidentParticleTrackInfo(
                                                CexmcTrackType  trackType ) :
    CexmcTrackInfo( trackType ), currentTrackLengthInTarget( 0. ),
    finalTrackLengthInTarget( 0. ), isStudiedProcessActivated( false ),
    needsTrackLengthResampling( false )
{
}


G4int  CexmcIncidentParticleTrackInfo::GetTypeInfo( void ) const
{
    return CexmcIncidentParticleTrackType;
}


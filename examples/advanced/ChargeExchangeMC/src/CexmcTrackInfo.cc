/*
 * ============================================================================
 *
 *       Filename:  CexmcTrackInfo.cc
 *
 *    Description:  track info
 *
 *        Version:  1.0
 *        Created:  22.11.2009 18:50:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcTrackInfo.hh"


CexmcTrackInfo::CexmcTrackInfo( CexmcTrackType  trackType, G4int  copyNumber ) :
    trackType( trackType ), copyNumber( copyNumber )
{
}


void  CexmcTrackInfo::Print( void ) const
{
}


G4int  CexmcTrackInfo::GetTypeInfo( void ) const
{
    return CexmcBasicTrackType;
}


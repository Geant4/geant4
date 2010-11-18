/*
 * ============================================================================
 *
 *       Filename:  CexmcEventInfo.cc
 *
 *    Description:  event information passed to run manager
 *
 *        Version:  1.0
 *        Created:  04.12.2009 15:50:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcEventInfo.hh"


CexmcEventInfo::CexmcEventInfo( G4bool  edTriggerIsOk, G4bool  tpTriggerIsOk,
                                G4bool  reconstructionIsOk ) :
    edTriggerIsOk( edTriggerIsOk ), tpTriggerIsOk( tpTriggerIsOk ),
    reconstructionIsOk( reconstructionIsOk )
{
}


void  CexmcEventInfo::Print( void ) const
{
}


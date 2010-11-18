/*
 * ============================================================================
 *
 *       Filename:  CexmcEventFastSObject.cc
 *
 *    Description:  event data serialization helper
 *
 *        Version:  1.0
 *        Created:  02.01.2010 20:31:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcEventFastSObject.hh"


CexmcEventFastSObject::CexmcEventFastSObject()
{
}


CexmcEventFastSObject::CexmcEventFastSObject( G4int  eventId,
                G4double  opCosThetaSCM, G4bool  edDigitizerHasTriggered,
                G4bool  edDigitizerMonitorHasTriggered ) :
    eventId( eventId ), opCosThetaSCM( opCosThetaSCM ),
    edDigitizerHasTriggered( edDigitizerHasTriggered ),
    edDigitizerMonitorHasTriggered( edDigitizerMonitorHasTriggered )
{
}


/*
 * ============================================================================
 *
 *       Filename:  CexmcEventSObject.cc
 *
 *    Description:  event data serialization helper
 *
 *        Version:  1.0
 *        Created:  30.12.2009 17:10:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcEventSObject.hh"
#include "CexmcTrackPointInfo.hh"


CexmcEventSObject::CexmcEventSObject() : edDigitizerMonitorHasTriggered( true )
{
}


CexmcEventSObject::CexmcEventSObject( G4int  eventId,
        G4bool  edDigitizerMonitorHasTriggered, G4double  monitorED,
        G4double  vetoCounterEDLeft, G4double  vetoCounterEDRight,
        G4double  calorimeterEDLeft, G4double  calorimeterEDRight,
        const CexmcEnergyDepositCalorimeterCollection &
                                                calorimeterEDLeftCollection,
        const CexmcEnergyDepositCalorimeterCollection &
                                                calorimeterEDRightCollection,
        const CexmcTrackPointInfo &  monitorTP,
        const CexmcTrackPointInfo &  targetTPBeamParticle,
        const CexmcTrackPointInfo &  targetTPOutputParticle,
        const CexmcTrackPointInfo &  targetTPNucleusParticle,
        const CexmcTrackPointInfo &
                                    targetTPOutputParticleDecayProductParticle1,
        const CexmcTrackPointInfo &
                                    targetTPOutputParticleDecayProductParticle2,
        const CexmcTrackPointInfo &  vetoCounterTPLeft,
        const CexmcTrackPointInfo &  vetoCounterTPRight,
        const CexmcTrackPointInfo &  calorimeterTPLeft,
        const CexmcTrackPointInfo &  calorimeterTPRight,
        const CexmcProductionModelData &  productionModelData ) :
    eventId( eventId ),
    edDigitizerMonitorHasTriggered( edDigitizerMonitorHasTriggered ),
    monitorED( monitorED ), vetoCounterEDLeft( vetoCounterEDLeft ),
    vetoCounterEDRight( vetoCounterEDRight ),
    calorimeterEDLeft( calorimeterEDLeft ),
    calorimeterEDRight( calorimeterEDRight ),
    calorimeterEDLeftCollection( calorimeterEDLeftCollection ),
    calorimeterEDRightCollection( calorimeterEDRightCollection ),
    monitorTP( monitorTP ), targetTPBeamParticle( targetTPBeamParticle ),
    targetTPOutputParticle( targetTPOutputParticle ),
    targetTPNucleusParticle( targetTPNucleusParticle ),
    targetTPOutputParticleDecayProductParticle1(
                                targetTPOutputParticleDecayProductParticle1 ),
    targetTPOutputParticleDecayProductParticle2(
                                targetTPOutputParticleDecayProductParticle2 ),
    vetoCounterTPLeft( vetoCounterTPLeft ),
    vetoCounterTPRight( vetoCounterTPRight ),
    calorimeterTPLeft( calorimeterTPLeft ),
    calorimeterTPRight( calorimeterTPRight ),
    productionModelData( productionModelData )
{
}


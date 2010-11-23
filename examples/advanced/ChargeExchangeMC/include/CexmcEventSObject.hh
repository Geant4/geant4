/*
 * =============================================================================
 *
 *       Filename:  CexmcEventSObject.hh
 *
 *    Description:  event data serialization helper
 *
 *        Version:  1.0
 *        Created:  30.12.2009 16:54:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_EVENT_SOBJECT_HH
#define CEXMC_EVENT_SOBJECT_HH

#ifdef CEXMC_USE_PERSISTENCY

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include "CexmcSimpleTrackPointInfoStore.hh"
#include "CexmcSimpleProductionModelDataStore.hh"
#include "CexmcCommon.hh"

class  CexmcTrackPointInfo;
class  CexmcProductionModelData;


class  CexmcEventSObject
{
    friend class  boost::serialization::access;
    friend class  CexmcRunManager;
#ifdef CEXMC_USE_CUSTOM_FILTER
    friend class  CexmcASTEval;
#endif

    public:
        CexmcEventSObject();

        CexmcEventSObject( G4int  eventId,
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
                const CexmcProductionModelData &  productionModelData );

    private:
        template  < typename  Archive >
        void  serialize( Archive &  archive, const unsigned int  version );

    private:
        G4int                                    eventId;

        G4bool                                   edDigitizerMonitorHasTriggered;

        G4double                                 monitorED;

        G4double                                 vetoCounterEDLeft;

        G4double                                 vetoCounterEDRight;

        G4double                                 calorimeterEDLeft;

        G4double                                 calorimeterEDRight;

        CexmcEnergyDepositCalorimeterCollection  calorimeterEDLeftCollection;

        CexmcEnergyDepositCalorimeterCollection  calorimeterEDRightCollection;

        CexmcSimpleTrackPointInfoStore           monitorTP;

        CexmcSimpleTrackPointInfoStore           targetTPBeamParticle;

        CexmcSimpleTrackPointInfoStore           targetTPOutputParticle;

        CexmcSimpleTrackPointInfoStore           targetTPNucleusParticle;

        CexmcSimpleTrackPointInfoStore
                                    targetTPOutputParticleDecayProductParticle1;

        CexmcSimpleTrackPointInfoStore
                                    targetTPOutputParticleDecayProductParticle2;

        CexmcSimpleTrackPointInfoStore           vetoCounterTPLeft;

        CexmcSimpleTrackPointInfoStore           vetoCounterTPRight;

        CexmcSimpleTrackPointInfoStore           calorimeterTPLeft;

        CexmcSimpleTrackPointInfoStore           calorimeterTPRight;

        CexmcSimpleProductionModelDataStore      productionModelData;
};


template  < typename  Archive >
void  CexmcEventSObject::serialize( Archive &  archive, const unsigned int )
{
    archive & eventId;
    archive & edDigitizerMonitorHasTriggered;
    archive & monitorED;
    archive & vetoCounterEDLeft;
    archive & vetoCounterEDRight;
    archive & calorimeterEDLeft;
    archive & calorimeterEDRight;
    archive & calorimeterEDLeftCollection;
    archive & calorimeterEDRightCollection;
    archive & monitorTP;
    archive & targetTPBeamParticle;
    archive & targetTPOutputParticle;
    archive & targetTPNucleusParticle;
    archive & targetTPOutputParticleDecayProductParticle1;
    archive & targetTPOutputParticleDecayProductParticle2;
    archive & vetoCounterTPLeft;
    archive & vetoCounterTPRight;
    archive & calorimeterTPLeft;
    archive & calorimeterTPRight;
    archive & productionModelData;
}

#endif

#endif


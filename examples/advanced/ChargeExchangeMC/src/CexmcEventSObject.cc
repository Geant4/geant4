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

#ifdef CEXMC_USE_PERSISTENCY

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

#endif


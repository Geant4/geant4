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

#include <boost/serialization/vector.hpp>
#include "CexmcSimpleTrackPointInfoStore.hh"
#include "CexmcSimpleProductionModelDataStore.hh"
#include "CexmcCommon.hh"


struct  CexmcEventSObject
{
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

    CexmcSimpleTrackPointInfoStore  targetTPOutputParticleDecayProductParticle1;

    CexmcSimpleTrackPointInfoStore  targetTPOutputParticleDecayProductParticle2;

    CexmcSimpleTrackPointInfoStore           vetoCounterTPLeft;

    CexmcSimpleTrackPointInfoStore           vetoCounterTPRight;

    CexmcSimpleTrackPointInfoStore           calorimeterTPLeft;

    CexmcSimpleTrackPointInfoStore           calorimeterTPRight;

    CexmcSimpleProductionModelDataStore      productionModelData;

    template  < typename  Archive >
    void  serialize( Archive &  archive, const unsigned int  version );
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


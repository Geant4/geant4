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
 *       Filename:  CexmcRunSObject.hh
 *
 *    Description:  run data serialization helper
 *
 *        Version:  1.0
 *        Created:  23.12.2009 14:24:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_RUN_SOBJECT_HH
#define CEXMC_RUN_SOBJECT_HH

#ifdef CEXMC_USE_PERSISTENCY

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include "CexmcSimpleDecayTableStore.hh"
#include "CexmcSimpleThreeVectorStore.hh"
#include "CexmcAngularRange.hh"
#include "CexmcSimpleRangeWithValue.hh"
#include "CexmcRun.hh"
#include "CexmcCommon.hh"


#define CEXMC_RUN_SOBJECT_VERSION 4


struct  CexmcRunSObject
{
    CexmcBasePhysicsUsed                 basePhysicsUsed;

    CexmcProductionModelType             productionModelType;

    std::string                          gdmlFileName;

    CexmcSimpleDecayTableStore           etaDecayTable;

    CexmcAngularRangeList                angularRanges;

    G4bool                               fermiMotionIsOn;

    std::vector< G4double >              calorimeterRegCuts;

    CexmcEventCountPolicy                eventCountPolicy;

    std::string                          beamParticle;

    CexmcSimpleThreeVectorStore          beamPos;

    CexmcSimpleThreeVectorStore          beamDir;

    G4double                             beamMomentumAmp;

    G4double                             beamFwhmPosX;

    G4double                             beamFwhmPosY;

    G4double                             beamFwhmDirX;

    G4double                             beamFwhmDirY;

    G4double                             beamFwhmMomentumAmp;

    G4double                             monitorEDThreshold;

    G4double                             vetoCounterEDLeftThreshold;

    G4double                             vetoCounterEDRightThreshold;

    G4double                             calorimeterEDLeftThreshold;

    G4double                             calorimeterEDRightThreshold;

    CexmcCalorimeterTriggerAlgorithm     calorimeterTriggerAlgorithm;

    CexmcOuterCrystalsVetoAlgorithm      outerCrystalsVetoAlgorithm;

    G4double                             outerCrystalsVetoFraction;

    G4bool                               applyFiniteCrystalResolution;

    CexmcEnergyRangeWithDoubleValueList  crystalResolutionData;

    CexmcCalorimeterEntryPointDefinitionAlgorithm
                                         epDefinitionAlgorithm;

    CexmcCalorimeterEntryPointDepthDefinitionAlgorithm
                                         epDepthDefinitionAlgorithm;

    CexmcCrystalSelectionAlgorithm       csAlgorithm;

    G4bool                               useInnerRefCrystal;

    G4double                             epDepth;

    G4bool                               useTableMass;

    G4bool                               useMassCut;

    G4double                             mCutOPCenter;

    G4double                             mCutNOPCenter;

    G4double                             mCutOPWidth;

    G4double                             mCutNOPWidth;

    G4double                             mCutAngle;

    G4bool                               useAbsorbedEnergyCut;

    G4double                             aeCutCLCenter;

    G4double                             aeCutCRCenter;

    G4double                             aeCutCLWidth;

    G4double                             aeCutCRWidth;

    G4double                             aeCutAngle;

    CexmcNmbOfHitsInRanges               nmbOfHitsSampled;

    CexmcNmbOfHitsInRanges               nmbOfHitsSampledFull;

    CexmcNmbOfHitsInRanges               nmbOfHitsTriggeredRealRange;

    CexmcNmbOfHitsInRanges               nmbOfHitsTriggeredRecRange;

    CexmcNmbOfHitsInRanges               nmbOfOrphanHits;

    G4int                                nmbOfFalseHitsTriggeredEDT;

    G4int                                nmbOfFalseHitsTriggeredRec;

    G4int                                nmbOfSavedEvents;

    G4int                                nmbOfSavedFastEvents;

    G4int                                numberOfEventsProcessed;

    G4int                                numberOfEventsProcessedEffective;

    G4int                                numberOfEventsToBeProcessed;

    std::string                          rProject;

    G4bool                               interactionsWithoutEDTWereSkipped;

    std::string                          cfFileName;

    CexmcEventDataVerboseLevel           evDataVerboseLevel;

    G4double                             proposedMaxIL;

    G4double                             expectedMomentumAmp;

    CexmcEDCollectionAlgoritm            edCollectionAlgorithm;

    unsigned int                         actualVersion;

    template  < typename  Archive >
    void  serialize( Archive &  archive, const unsigned int  version );
};


template  < typename  Archive >
void  CexmcRunSObject::serialize( Archive &  archive,
                                  const unsigned int  version )
{
    archive & basePhysicsUsed;
    archive & productionModelType;
    archive & gdmlFileName;
    archive & etaDecayTable;
    archive & angularRanges;
    archive & fermiMotionIsOn;
    archive & calorimeterRegCuts;
    archive & eventCountPolicy;
    archive & beamParticle;
    archive & beamPos;
    archive & beamDir;
    archive & beamMomentumAmp;
    archive & beamFwhmPosX;
    archive & beamFwhmPosY;
    archive & beamFwhmDirX;
    archive & beamFwhmDirY;
    archive & beamFwhmMomentumAmp;
    archive & monitorEDThreshold;
    archive & vetoCounterEDLeftThreshold;
    archive & vetoCounterEDRightThreshold;
    archive & calorimeterEDLeftThreshold;
    archive & calorimeterEDRightThreshold;
    archive & calorimeterTriggerAlgorithm;
    archive & outerCrystalsVetoAlgorithm;
    archive & outerCrystalsVetoFraction;
    archive & applyFiniteCrystalResolution;
    archive & crystalResolutionData;
    archive & epDefinitionAlgorithm;
    archive & epDepthDefinitionAlgorithm;
    archive & csAlgorithm;
    if ( version > 0 )
        archive & useInnerRefCrystal;
    archive & epDepth;
    archive & useTableMass;
    archive & useMassCut;
    archive & mCutOPCenter;
    archive & mCutNOPCenter;
    archive & mCutOPWidth;
    archive & mCutNOPWidth;
    archive & mCutAngle;
    archive & useAbsorbedEnergyCut;
    archive & aeCutCLCenter;
    archive & aeCutCRCenter;
    archive & aeCutCLWidth;
    archive & aeCutCRWidth;
    archive & aeCutAngle;
    archive & nmbOfHitsSampled;
    archive & nmbOfHitsSampledFull;
    archive & nmbOfHitsTriggeredRealRange;
    archive & nmbOfHitsTriggeredRecRange;
    archive & nmbOfOrphanHits;
    archive & nmbOfFalseHitsTriggeredEDT;
    archive & nmbOfFalseHitsTriggeredRec;
    archive & nmbOfSavedEvents;
    archive & nmbOfSavedFastEvents;
    archive & numberOfEventsProcessed;
    archive & numberOfEventsProcessedEffective;
    archive & numberOfEventsToBeProcessed;
    if ( version > 1 )
    {
        archive & rProject;
        archive & interactionsWithoutEDTWereSkipped;
        archive & cfFileName;
        archive & evDataVerboseLevel;
    }
    if ( version > 2 )
        archive & proposedMaxIL;
    if ( version > 3 )
    {
        archive & expectedMomentumAmp;
        archive & edCollectionAlgorithm;
    }

    actualVersion = version;
}


BOOST_CLASS_VERSION( CexmcRunSObject, CEXMC_RUN_SOBJECT_VERSION )

#endif

#endif


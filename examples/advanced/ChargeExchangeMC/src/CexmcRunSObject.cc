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
 *       Filename:  CexmcRunSObject.cc
 *
 *    Description:  run data serialization helper
 *
 *        Version:  1.0
 *        Created:  23.12.2009 14:41:49
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifdef CEXMC_USE_PERSISTENCY

#include "CexmcRunSObject.hh"

CexmcRunSObject::CexmcRunSObject()
{
}


CexmcRunSObject::CexmcRunSObject(
        CexmcBasePhysicsUsed  basePhysicsUsed,
        CexmcProductionModelType  productionModelType,
        const std::string &  gdmlFileName,
        const CexmcSimpleDecayTableStore &  etaDecayTable,
        const CexmcAngularRangeList &  angularRanges, G4bool  fermiMotionIsOn,
        const std::vector< G4double > &  calorimeterRegCuts,
        CexmcEventCountPolicy  eventCountPolicy,
        const std::string &  beamParticle,
        const CexmcSimpleThreeVectorStore &  beamPos,
        const CexmcSimpleThreeVectorStore &  beamDir, G4double  beamMomentumAmp,
        G4double  beamFwhmPosX, G4double  beamFwhmPosY, G4double  beamFwhmDirX,
        G4double  beamFwhmDirY, G4double  beamFwhmMomentumAmp,
        G4double  monitorEDThreshold, G4double  vetoCounterEDLeftThreshold,
        G4double  vetoCounterEDRightThreshold,
        G4double  calorimeterEDLeftThreshold,
        G4double  calorimeterEDRightThreshold,
        CexmcCalorimeterTriggerAlgorithm  calorimeterTriggerAlgorithm,
        CexmcOuterCrystalsVetoAlgorithm  outerCrystalsVetoAlgorithm,
        G4double  outerCrystalsVetoFraction,
        G4bool  applyFiniteCrystalResolution,
        const CexmcEnergyRangeWithDoubleValueList &  crystalResolutionData,
        CexmcCalorimeterEntryPointDefinitionAlgorithm  epDefinitionAlgorithm,
        CexmcCalorimeterEntryPointDepthDefinitionAlgorithm
                                                    epDepthDefinitionAlgorithm,
        CexmcCrystalSelectionAlgorithm  csAlgorithm, G4bool  useInnerRefCrystal,
        G4double  epDepth, G4bool  useTableMass, G4bool  useMassCut,
        G4double  mCutOPCenter, G4double  mCutNOPCenter, G4double  mCutOPWidth,
        G4double  mCutNOPWidth, G4double  mCutAngle,
        G4bool  useAbsorbedEnergyCut, G4double  aeCutCLCenter,
        G4double  aeCutCRCenter, G4double  aeCutCLWidth, G4double  aeCutCRWidth,
        G4double  aeCutAngle, CexmcNmbOfHitsInRanges  nmbOfHitsSampled,
        CexmcNmbOfHitsInRanges  nmbOfHitsSampledFull,
        CexmcNmbOfHitsInRanges  nmbOfHitsTriggeredRealRange,
        CexmcNmbOfHitsInRanges  nmbOfHitsTriggeredRecRange,
        CexmcNmbOfHitsInRanges  nmbOfOrphanHits,
        G4int  nmbOfFalseHitsTriggeredEDT, G4int  nmbOfFalseHitsTriggeredRec,
        G4int  nmbOfSavedEvents, G4int  nmbOfSavedFastEvents,
        G4int  numberOfEventsProcessed, G4int  numberOfEventsProcessedEffective,
        G4int  numberOfEventsToBeProcessed, const std::string &  rProject,
        G4bool  interactionsWithoutEDTWereSkipped,
        const std::string &  cfFileName,
        CexmcEventDataVerboseLevel  evDataVerboseLevel,
        G4double  proposedMaxIL ) :
    basePhysicsUsed( basePhysicsUsed ),
    productionModelType( productionModelType ), gdmlFileName( gdmlFileName ),
    etaDecayTable( etaDecayTable ), angularRanges( angularRanges ),
    fermiMotionIsOn( fermiMotionIsOn ),
    calorimeterRegCuts( calorimeterRegCuts ),
    eventCountPolicy( eventCountPolicy ), beamParticle( beamParticle ),
    beamPos( beamPos ), beamDir( beamDir ), beamMomentumAmp( beamMomentumAmp ),
    beamFwhmPosX( beamFwhmPosX ), beamFwhmPosY( beamFwhmPosY ),
    beamFwhmDirX( beamFwhmDirX ), beamFwhmDirY( beamFwhmDirY ),
    beamFwhmMomentumAmp( beamFwhmMomentumAmp ),
    monitorEDThreshold( monitorEDThreshold ),
    vetoCounterEDLeftThreshold( vetoCounterEDLeftThreshold ),
    vetoCounterEDRightThreshold( vetoCounterEDRightThreshold ),
    calorimeterEDLeftThreshold( calorimeterEDLeftThreshold ),
    calorimeterEDRightThreshold( calorimeterEDRightThreshold ),
    calorimeterTriggerAlgorithm( calorimeterTriggerAlgorithm ),
    outerCrystalsVetoAlgorithm( outerCrystalsVetoAlgorithm ),
    outerCrystalsVetoFraction( outerCrystalsVetoFraction ),
    applyFiniteCrystalResolution( applyFiniteCrystalResolution ),
    crystalResolutionData( crystalResolutionData ),
    epDefinitionAlgorithm( epDefinitionAlgorithm ),
    epDepthDefinitionAlgorithm( epDepthDefinitionAlgorithm ),
    csAlgorithm( csAlgorithm ), useInnerRefCrystal( useInnerRefCrystal ),
    epDepth( epDepth ), useTableMass( useTableMass ), useMassCut( useMassCut ),
    mCutOPCenter( mCutOPCenter ), mCutNOPCenter( mCutNOPCenter ),
    mCutOPWidth( mCutOPWidth ), mCutNOPWidth( mCutNOPWidth ),
    mCutAngle( mCutAngle ), useAbsorbedEnergyCut( useAbsorbedEnergyCut ),
    aeCutCLCenter( aeCutCLCenter ), aeCutCRCenter( aeCutCRCenter ),
    aeCutCLWidth( aeCutCLWidth ), aeCutCRWidth( aeCutCRWidth ),
    aeCutAngle( aeCutAngle ), nmbOfHitsSampled( nmbOfHitsSampled ),
    nmbOfHitsSampledFull( nmbOfHitsSampledFull ),
    nmbOfHitsTriggeredRealRange( nmbOfHitsTriggeredRealRange ),
    nmbOfHitsTriggeredRecRange( nmbOfHitsTriggeredRecRange ),
    nmbOfOrphanHits( nmbOfOrphanHits ),
    nmbOfFalseHitsTriggeredEDT( nmbOfFalseHitsTriggeredEDT ),
    nmbOfFalseHitsTriggeredRec( nmbOfFalseHitsTriggeredRec ),
    nmbOfSavedEvents( nmbOfSavedEvents ),
    nmbOfSavedFastEvents( nmbOfSavedFastEvents ),
    numberOfEventsProcessed( numberOfEventsProcessed ),
    numberOfEventsProcessedEffective( numberOfEventsProcessedEffective ),
    numberOfEventsToBeProcessed( numberOfEventsToBeProcessed ),
    rProject( rProject ),
    interactionsWithoutEDTWereSkipped( interactionsWithoutEDTWereSkipped ),
    cfFileName( cfFileName ), evDataVerboseLevel( evDataVerboseLevel ),
    proposedMaxIL( proposedMaxIL ), actualVersion( CEXMC_RUN_SOBJECT_VERSION )
{
}

#endif


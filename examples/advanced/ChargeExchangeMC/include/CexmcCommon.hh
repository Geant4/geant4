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
 *       Filename:  CexmcCommon.hh
 *
 *    Description:  common declarations
 *
 *        Version:  1.0
 *        Created:  01.11.2009 00:09:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_COMMON_HH
#define CEXMC_COMMON_HH

#include <vector>
#include <limits>
#include <G4String.hh>
#include <G4Types.hh>

#define CEXMC_LINE_START  "--- Cexmc ---  "


typedef std::vector< G4double >  CexmcEnergyDepositCrystalRowCollection;

typedef std::vector< CexmcEnergyDepositCrystalRowCollection >
                                 CexmcEnergyDepositCalorimeterCollection;


const G4double  CexmcDblMax( std::numeric_limits< double >::max() );

const G4String  CexmcStudiedProcessFirstName( "studiedProcess_" );
const G4String  CexmcStudiedProcessLastName( "Cexmc" );
const G4String  CexmcStudiedProcessFullName( CexmcStudiedProcessFirstName +
                                             CexmcStudiedProcessLastName );
const G4String  CexmcChargeExchangeProductionModelName( "ChargeExchange" );
const G4String  CexmcChargeExchangeInteractionName( "Cexmc" +
                                    CexmcChargeExchangeProductionModelName );

const G4String  CexmcEDDigitizerName( "EDDig" );
const G4String  CexmcTPDigitizerName( "TPDig" );

const G4String  CexmcScenePrimitivesDescription( "CexmcScenePrimitives" );

const G4double  CexmcFwhmToStddev( 0.42466 );
const G4double  CexmcInvalidCosTheta( 2.0 );

const G4int     CexmcInvalidTrackId( -1 );


enum CexmcBasePhysicsUsed
{
    CexmcNoBasePhysics,
    Cexmc_QGSP_BERT,
    Cexmc_QGSP_BIC_EMY,
    Cexmc_FTFP_BERT
};


enum  CexmcProductionModelType
{
    CexmcUnknownProductionModel,
    CexmcPionZeroProduction,
    CexmcEtaProduction
};


enum  CexmcTriggerType
{
    CexmcTPT,
    CexmcEDT,
    CexmcRT
};


enum  CexmcEventCountPolicy
{
    CexmcCountAllEvents,
    CexmcCountEventsWithInteraction,
    CexmcCountEventsWithTrigger
};


enum  CexmcTrackType
{
    CexmcInsipidTrack,
    CexmcBeamParticleTrack,
    CexmcOutputParticleTrack,
    CexmcNucleusParticleTrack,
    CexmcOutputParticleDecayProductTrack
};


enum  CexmcTrackTypeInfo
{
    CexmcBasicTrackType,
    CexmcIncidentParticleTrackType
};


enum  CexmcSide
{
    CexmcLeft,
    CexmcRight
};


enum  CexmcOuterCrystalsVetoAlgorithm
{
    CexmcNoOuterCrystalsVeto,
    CexmcMaximumEDInASingleOuterCrystalVeto,
    CexmcFractionOfEDInOuterCrystalsVeto
};


enum  CexmcCalorimeterTriggerAlgorithm
{
    CexmcAllCrystalsMakeEDTriggerThreshold,
    CexmcInnerCrystalsMakeEDTriggerThreshold
};


enum  CexmcCalorimeterEntryPointDefinitionAlgorithm
{
    CexmcEntryPointInTheCenter,
    CexmcEntryPointInTheCenterOfCrystalWithMaxED,
    CexmcEntryPointByLinearEDWeights,
    CexmcEntryPointBySqrtEDWeights
};


enum  CexmcCalorimeterEntryPointDepthDefinitionAlgorithm
{
    CexmcEntryPointDepthPlain,
    CexmcEntryPointDepthSphere
};


enum  CexmcCrystalSelectionAlgorithm
{
    CexmcSelectAllCrystals,
    CexmcSelectAdjacentCrystals
};


enum  CexmcEDCollectionAlgoritm
{
    CexmcCollectEDInAllCrystals,
    CexmcCollectEDInAdjacentCrystals
};


enum  CexmcEventDataVerboseLevel
{
    CexmcWriteNoEventData,
    CexmcWriteEventDataOnEveryEDT,
    CexmcWriteEventDataOnEveryTPT
};


enum  CexmcOutputDataType
{
    CexmcOutputRun,
    CexmcOutputGeometry,
    CexmcOutputEvents
};


#endif


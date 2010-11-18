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

const G4double  CexmcFwhmToStddev( 0.42466 );
const G4double  CexmcInvalidCosTheta( 2.0 );

const G4int     CexmcInvalidTrackId( -1 );


enum CexmcBasePhysicsUsed
{
    CexmcNoBasePhysics,
    Cexmc_QGSP_BERT,
    Cexmc_QGSP_BIC_EMY
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


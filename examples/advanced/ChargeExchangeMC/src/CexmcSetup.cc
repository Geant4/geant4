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
 *       Filename:  CexmcSetup.cc
 *
 *    Description:  physical setup
 *
 *        Version:  1.0
 *        Created:  10.10.2009 23:00:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <G4GDMLParser.hh>
#include <G4MultiFunctionalDetector.hh>
#include <G4SDManager.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4Box.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4ProductionCuts.hh>
#include <G4VUserPhysicsList.hh>
#include <G4SystemOfUnits.hh>
#include "CexmcSetup.hh"
#include "CexmcPrimitiveScorer.hh"
#include "CexmcTrackPoints.hh"
#include "CexmcTrackPointsInLeftRightSet.hh"
#include "CexmcTrackPointsInCalorimeter.hh"
#include "CexmcTrackPointsFilter.hh"
#include "CexmcSimpleEnergyDeposit.hh"
#include "CexmcEnergyDepositInLeftRightSet.hh"
#include "CexmcEnergyDepositInCalorimeter.hh"
#include "CexmcRunManager.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcException.hh"


CexmcSetup::CexmcSetup( const G4String &  gdmlFile_,
                        G4bool  validateGDMLFile_ ) :
    world( 0 ), gdmlFile( gdmlFile_ ), validateGDMLFile( validateGDMLFile_ ),
    calorimeterRegionInitialized( false ),
    calorimeterGeometryDataInitialized( false ), monitorVolume( NULL ),
    vetoCounterVolume( NULL ), calorimeterVolume( NULL ), targetVolume( NULL ),
    rightVetoCounter( NULL ), rightCalorimeter( NULL )
{
}


G4VPhysicalVolume *  CexmcSetup::Construct( void )
{
    if ( world )
        return world;

    G4GDMLParser  gdmlParser;

    gdmlParser.Read( gdmlFile, validateGDMLFile );
    world = gdmlParser.GetWorldVolume();

    SetupSpecialVolumes( gdmlParser );

    ReadTransforms( gdmlParser );

    ReadRightDetectors();

    CexmcRunManager *  runManager( static_cast< CexmcRunManager * >(
                                            G4RunManager::GetRunManager() ) );

    runManager->SetupConstructionHook();

    const CexmcPhysicsManager *  physicsManager(
            dynamic_cast< const CexmcPhysicsManager * >(
                                        runManager->GetUserPhysicsList() ) );

    if ( ! physicsManager )
        throw CexmcException( CexmcWeirdException );

    CexmcPhysicsManager *        thePhysicsManager(
            const_cast< CexmcPhysicsManager * >( physicsManager ) );
    thePhysicsManager->SetupConstructionHook( this );

    return world;
}


void  CexmcSetup::SetupSpecialVolumes( const G4GDMLParser &  gdmlParser )
{
    G4MultiFunctionalDetector *   detector[ CexmcNumberOfDetectorRoles ] =
                                                                    { NULL };
    const G4LogicalVolumeStore *  lvs( G4LogicalVolumeStore::GetInstance() );

    for ( std::vector< G4LogicalVolume * >::const_iterator
                        lvIter( lvs->begin() ); lvIter != lvs->end(); ++lvIter )
    {
        G4String           volumeName( G4String( ( *lvIter )->GetName() ) );
        G4GDMLAuxListType  auxInfo( gdmlParser.GetVolumeAuxiliaryInformation(
                                                                    *lvIter ) );
        CexmcDetectorRole  curDetectorRole( CexmcNumberOfDetectorRoles );

        for ( G4GDMLAuxListType::const_iterator  pair( auxInfo.begin() );
                                              pair != auxInfo.end(); ++pair )
        {
            CexmcPrimitiveScorer *  scorer( NULL );
            G4String                detectorName( "uninitialized" );
            do
            {
                if ( pair->type == "EnergyDepositDetector" )
                {
                    do
                    {
                        if ( pair->value == "MonitorRole" )
                        {
                            AssertAndAsignDetectorRole( curDetectorRole,
                                                CexmcMonitorDetectorRole );
                            scorer = new CexmcSimpleEnergyDeposit(
                                    CexmcDetectorTypeName[ CexmcEDDetector ] );
                            break;
                        }
                        if ( pair->value == "VetoCounterRole" )
                        {
                            AssertAndAsignDetectorRole( curDetectorRole,
                                                CexmcVetoCounterDetectorRole );
                            scorer = new CexmcEnergyDepositInLeftRightSet(
                                    CexmcDetectorTypeName[ CexmcEDDetector ],
                                    this );
                            break;
                        }
                        if ( pair->value == "CalorimeterRole" )
                        {
                            AssertAndAsignDetectorRole( curDetectorRole,
                                                CexmcCalorimeterDetectorRole );
                            scorer = new CexmcEnergyDepositInCalorimeter(
                                    CexmcDetectorTypeName[ CexmcEDDetector ],
                                    this );
                            break;
                        }
                    } while ( false );
                    detectorName = CexmcDetectorRoleName[ curDetectorRole ];
                    G4cout << CEXMC_LINE_START "ED Scorer of detector role '" <<
                               detectorName << "' in volume '" << volumeName <<
                               "'" << G4endl;
                    break;
                }
                if ( pair->type == "TrackPointsDetector" )
                {
                    do
                    {
                        if ( pair->value == "MonitorRole" )
                        {
                            AssertAndAsignDetectorRole( curDetectorRole,
                                                CexmcMonitorDetectorRole );
                            scorer = new CexmcTrackPoints(
                                    CexmcDetectorTypeName[ CexmcTPDetector ] );
                            break;
                        }
                        if ( pair->value == "VetoCounterRole"  )
                        {
                            AssertAndAsignDetectorRole( curDetectorRole,
                                                CexmcVetoCounterDetectorRole );
                            scorer = new CexmcTrackPointsInLeftRightSet(
                                    CexmcDetectorTypeName[ CexmcTPDetector ],
                                    this );
                            break;
                        }
                        if ( pair->value == "CalorimeterRole" )
                        {
                            AssertAndAsignDetectorRole( curDetectorRole,
                                                CexmcCalorimeterDetectorRole );
                            scorer = new CexmcTrackPointsInCalorimeter(
                                    CexmcDetectorTypeName[ CexmcTPDetector ],
                                    this );
                            break;
                        }
                        if ( pair->value == "TargetRole" )
                        {
                            AssertAndAsignDetectorRole( curDetectorRole,
                                                CexmcTargetDetectorRole );
                            scorer = new CexmcTrackPoints(
                                    CexmcDetectorTypeName[ CexmcTPDetector ] );
                            break;
                        }
                    } while ( false );
                    detectorName = CexmcDetectorRoleName[ curDetectorRole ];
                    G4cout << CEXMC_LINE_START "TP Scorer of detector role '" <<
                               detectorName << "' in volume '" << volumeName <<
                               "'" << G4endl;
                    if ( scorer )
                    {
                        CexmcTrackPointsFilter *  filter(
                                 new CexmcTrackPointsFilter( "trackPoints" ) );
                        scorer->SetFilter( filter );
                    }
                    break;
                }
                if ( pair->type == "SensitiveRegion" )
                {
                    do
                    {
                        if ( pair->value == "CalorimeterRegion" )
                        {
                            G4Region *  region( NULL );
                            if ( calorimeterRegionInitialized )
                            {
                                region = G4RegionStore::GetInstance()->
                                        GetRegion( CexmcCalorimeterRegionName );
                            }
                            else
                            {
                                region = new G4Region(
                                                CexmcCalorimeterRegionName );
                                G4ProductionCuts *  cuts(
                                                        new G4ProductionCuts );
                                G4double  defaultProductionCut( 1.0 * mm );
                                const G4VUserPhysicsList *  physicsList(
                                            G4RunManager::GetRunManager()->
                                                        GetUserPhysicsList() );
                                if ( physicsList )
                                    defaultProductionCut =
                                            physicsList->GetDefaultCutValue();
                                cuts->SetProductionCut( defaultProductionCut );
                                region->SetProductionCuts( cuts );
                                calorimeterRegionInitialized = true;
                            }
                            region->AddRootLogicalVolume( *lvIter );
                            break;
                        }
                    } while ( false );
                    G4cout << CEXMC_LINE_START "Sensitive Region for logical "
                               "volume '" << volumeName << "' registered" <<
                               G4endl;
                    break;
                }
                if ( pair->type == "SpecialVolume" )
                {
                    do
                    {
                        if ( pair->value == "Monitor" )
                        {
                            monitorVolume = *lvIter;
                            G4cout << CEXMC_LINE_START "Monitor volume '";
                            break;
                        }
                        if ( pair->value == "VetoCounter" )
                        {
                            vetoCounterVolume = *lvIter;
                            G4cout << CEXMC_LINE_START "VetoCounter volume '";
                            break;
                        }
                        if ( pair->value == "Calorimeter" )
                        {
                            calorimeterVolume = *lvIter;
                            G4cout << CEXMC_LINE_START "Calorimeter volume '";
                            ReadCalorimeterGeometryData( *lvIter );
                            calorimeterGeometryDataInitialized = true;
                            break;
                        }
                        if ( pair->value == "Target" )
                        {
                            targetVolume = *lvIter;
                            G4cout << CEXMC_LINE_START "Target volume '";
                            break;
                        }
                    } while ( false );
                    G4cout << volumeName << "' registered" << G4endl;
                    break;
                }
            }
            while ( false );

            if ( scorer )
            {
                /* curDetectorRole must be intact when scorer is not NULL */
                if ( ! detector[ curDetectorRole ] )
                {
                    detector[ curDetectorRole ] =
                                new G4MultiFunctionalDetector( detectorName );
                }
                detector[ curDetectorRole ]->RegisterPrimitive( scorer );
                /* now that scorer has initialized pointer to its detector, its
                 * messenger's path shall be properly initialized as well */
                scorer->InitializeMessenger();
                /* NB: logical volumes in GDML file may not have multiple
                 * detector roles: for example volume Monitor may have only one
                 * role MonitorRole. This restriction arises from that fact that
                 * a logical volume may contain only one sensitive detector. */
                ( *lvIter )->SetSensitiveDetector(
                                                detector[ curDetectorRole ] );
            }
        }
    }

    if ( ! calorimeterRegionInitialized )
        throw CexmcException( CexmcCalorimeterRegionNotInitialized );

    if ( ! calorimeterGeometryDataInitialized )
        throw CexmcException( CexmcCalorimeterGeometryDataNotInitialized );

    for ( G4int  i( 0 ); i < CexmcNumberOfDetectorRoles; ++i )
    {
        if ( detector[ i ] )
            G4SDManager::GetSDMpointer()->AddNewDetector( detector[ i ] );
    }
}


void  CexmcSetup::ReadTransforms( const G4GDMLParser &  gdmlParser )
{
    G4ThreeVector     position( gdmlParser.GetPosition( "TargetPos" ) );
    G4ThreeVector     rotation( gdmlParser.GetRotation( "TargetRot" ) );
    G4RotationMatrix  rm;

    RotateMatrix( rotation, rm );
    targetTransform.SetNetTranslation( position );
    targetTransform.SetNetRotation( rm );

    position = gdmlParser.GetPosition( "CalorimeterLeftPos" );
    rotation = gdmlParser.GetRotation( "CalorimeterLeftRot" );
    rm = G4RotationMatrix();
    RotateMatrix( rotation, rm );
    calorimeterLeftTransform.SetNetTranslation( position );
    calorimeterLeftTransform.SetNetRotation( rm );

    position = gdmlParser.GetPosition( "CalorimeterRightPos" );
    rotation = gdmlParser.GetRotation( "CalorimeterRightRot" );
    rm = G4RotationMatrix();
    RotateMatrix( rotation, rm );
    calorimeterRightTransform.SetNetTranslation( position );
    calorimeterRightTransform.SetNetRotation( rm );
}


void  CexmcSetup::ReadCalorimeterGeometryData(
                                            const G4LogicalVolume *  lVolume )
{
    if ( lVolume->GetNoDaughters() == 0 )
        throw CexmcException( CexmcIncompatibleGeometry );

    G4VPhysicalVolume *  pVolume( lVolume->GetDaughter( 0 ) );
    EAxis                axis;
    G4double             width;
    G4double             offset;
    G4bool               consuming;

    if ( ! pVolume )
        throw CexmcException( CexmcIncompatibleGeometry );

    if ( pVolume->IsReplicated() )
    {
        pVolume->GetReplicationData( axis,
                                     calorimeterGeometry.nCrystalsInColumn,
                                     width, offset, consuming );
    }

    lVolume = pVolume->GetLogicalVolume();

    if ( lVolume->GetNoDaughters() == 0 )
        throw CexmcException( CexmcIncompatibleGeometry );

    pVolume = lVolume->GetDaughter( 0 );

    if ( ! pVolume )
        throw CexmcException( CexmcIncompatibleGeometry );

    if ( pVolume->IsReplicated() )
    {
        pVolume->GetReplicationData( axis, calorimeterGeometry.nCrystalsInRow,
                                     width, offset, consuming );
    }

    lVolume = pVolume->GetLogicalVolume();

    /* NB: this is not necessarily a crystal itself as far as crystals can be
     * wrapped in paper and other materials, but this is what reconstructor and
     * digitizers really need */
    G4Box *  crystalBox( dynamic_cast< G4Box * >( lVolume->GetSolid() ) );

    if ( ! crystalBox )
        throw CexmcException( CexmcIncompatibleGeometry );

    calorimeterGeometry.crystalWidth = crystalBox->GetXHalfLength() * 2;
    calorimeterGeometry.crystalHeight = crystalBox->GetYHalfLength() * 2;
    calorimeterGeometry.crystalLength = crystalBox->GetZHalfLength() * 2;
}


void  CexmcSetup::ConvertToCrystalGeometry( const G4ThreeVector &  src,
                    G4int &  row, G4int &  column, G4ThreeVector &  dst ) const
{
    G4int     nCrystalsInColumn( calorimeterGeometry.nCrystalsInColumn );
    G4int     nCrystalsInRow( calorimeterGeometry.nCrystalsInRow );
    G4double  crystalWidth( calorimeterGeometry.crystalWidth );
    G4double  crystalHeight( calorimeterGeometry.crystalHeight );

    row = G4int( ( src.y() + crystalHeight * nCrystalsInColumn / 2 ) /
                 crystalHeight );
    column = G4int( ( src.x() + crystalWidth * nCrystalsInRow / 2 ) /
                    crystalWidth );
    G4double   xInCalorimeterOffset(
                    ( G4double( column ) - G4double( nCrystalsInRow ) / 2 ) *
                                            crystalWidth + crystalWidth / 2 );
    G4double   yInCalorimeterOffset(
                    ( G4double( row ) - G4double( nCrystalsInColumn ) / 2 ) *
                                            crystalHeight + crystalHeight / 2 );
    dst.setX( src.x() - xInCalorimeterOffset );
    dst.setY( src.y() - yInCalorimeterOffset );
}


void  CexmcSetup::ReadRightDetectors( void )
{
    G4PhysicalVolumeStore *  pvs( G4PhysicalVolumeStore::GetInstance() );

    for ( std::vector< G4VPhysicalVolume * >::const_iterator  k( pvs->begin() );
                                                        k != pvs->end(); ++k )
    {
        /* FIXME: it would be more reasonable to find detectors from volumes
         * tagged with 'EnergyDepositDetector' or 'TrackPointsDetector', and not
         * from volumes tagged with 'SpecialVolume' as it is done here. However
         * in case of calorimeters the role of detectors are played by crystal
         * volumes, not the calorimeters themselves, but only calorimeters can
         * hold information about their left or right positions! Thus, following
         * considerations of convenience, right detectors for all detector roles
         * are chosen from volumes tagged with 'SpecialVolume' */
        do
        {
            if ( ( *k )->GetLogicalVolume() == vetoCounterVolume )
            {
                if ( ( *k )->GetName().contains( "Right" ) )
                    rightVetoCounter = *k;
                break;
            }
            if ( ( *k )->GetLogicalVolume() == calorimeterVolume )
            {
                if ( ( *k )->GetName().contains( "Right" ) )
                    rightCalorimeter = *k;
                break;
            }
        } while ( false );
    }
}


void  CexmcSetup::AssertAndAsignDetectorRole( CexmcDetectorRole &  detectorRole,
                                              CexmcDetectorRole  value )
{
    if ( detectorRole != CexmcNumberOfDetectorRoles && detectorRole != value )
        throw CexmcException( CexmcMultipleDetectorRoles );

    detectorRole = value;
}


void  CexmcSetup::RotateMatrix( const G4ThreeVector &  rot,
                                G4RotationMatrix &  rm )
{
    rm.rotateX( rot.x() );
    rm.rotateY( rot.y() );
    rm.rotateZ( rot.z() );
}


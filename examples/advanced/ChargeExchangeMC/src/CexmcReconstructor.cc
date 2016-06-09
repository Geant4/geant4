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
 *       Filename:  CexmcReconstructor.cc
 *
 *    Description:  reconstructor base class
 *
 *        Version:  1.0
 *        Created:  02.12.2009 16:21:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcReconstructor.hh"
#include "CexmcReconstructorMessenger.hh"
#include "CexmcEnergyDepositStore.hh"
#include "CexmcRunManager.hh"


CexmcReconstructor::CexmcReconstructor() : hasBasicTrigger( false ),
    epDefinitionAlgorithm( CexmcEntryPointBySqrtEDWeights ),
    epDepthDefinitionAlgorithm( CexmcEntryPointDepthPlain ),
    csAlgorithm( CexmcSelectAllCrystals ), useInnerRefCrystal( false ),
    epDepth( 0 ), theAngle( 0 ), calorimeterEDLeftAdjacent( 0 ),
    calorimeterEDRightAdjacent( 0 ), collectEDInAdjacentCrystals( false ),
    targetEPInitialized( false ), messenger( NULL )
{
    G4RunManager *      runManager( G4RunManager::GetRunManager() );
    const CexmcSetup *  setup( static_cast< const CexmcSetup * >(
                                runManager->GetUserDetectorConstruction() ) );
    calorimeterGeometry = setup->GetCalorimeterGeometry();
    targetTransform = setup->GetTargetTransform();
    calorimeterLeftTransform = setup->GetCalorimeterLeftTransform();
    calorimeterRightTransform = setup->GetCalorimeterRightTransform();

    messenger = new CexmcReconstructorMessenger( this );
}


CexmcReconstructor::~CexmcReconstructor()
{
    delete messenger;
}


void  CexmcReconstructor::Reconstruct(
                                    const CexmcEnergyDepositStore *  edStore )
{
    ReconstructEntryPoints( edStore );
    if ( hasBasicTrigger )
        ReconstructTargetPoint();
    if ( hasBasicTrigger )
        ReconstructAngle();
}


G4bool  CexmcReconstructor::HasFullTrigger( void ) const
{
    return hasBasicTrigger;
}


void  CexmcReconstructor::ReconstructEntryPoints(
                                    const CexmcEnergyDepositStore *  edStore )
{
    G4int     columnLeft( edStore->calorimeterEDLeftMaxX );
    G4int     rowLeft( edStore->calorimeterEDLeftMaxY );
    G4int     columnRight( edStore->calorimeterEDRightMaxX );
    G4int     rowRight( edStore->calorimeterEDRightMaxY );
    G4double  crystalLength( calorimeterGeometry.crystalLength );

    if ( useInnerRefCrystal )
    {
        TransformToAdjacentInnerCrystal( columnLeft, rowLeft );
        TransformToAdjacentInnerCrystal( columnRight, rowRight );
    }

    calorimeterEPLeftPosition.setX( 0 );
    calorimeterEPLeftPosition.setY( 0 );
    calorimeterEPLeftPosition.setZ( -crystalLength / 2 + epDepth );
    calorimeterEPLeftDirection.setX( 0 );
    calorimeterEPLeftDirection.setY( 0 );
    calorimeterEPLeftDirection.setZ( 0 );
    calorimeterEPRightPosition.setX( 0 );
    calorimeterEPRightPosition.setY( 0 );
    calorimeterEPRightPosition.setZ( -crystalLength / 2 + epDepth );
    calorimeterEPRightDirection.setX( 0 );
    calorimeterEPRightDirection.setY( 0 );
    calorimeterEPRightDirection.setZ( 0 );

    G4bool  edInAdjacentCrystalsCollected( false );

    switch ( epDefinitionAlgorithm )
    {
    case CexmcEntryPointInTheCenter :
        break;
    case CexmcEntryPointInTheCenterOfCrystalWithMaxED :
        {
            G4int     nCrystalsInColumn(
                                    calorimeterGeometry.nCrystalsInColumn );
            G4int     nCrystalsInRow( calorimeterGeometry.nCrystalsInRow );
            G4double  crystalWidth( calorimeterGeometry.crystalWidth );
            G4double  crystalHeight( calorimeterGeometry.crystalHeight );

            calorimeterEPLeftPosition.setX(
                ( G4double( columnLeft ) - G4double( nCrystalsInRow ) / 2 ) *
                                        crystalWidth  + crystalWidth / 2 );
            calorimeterEPLeftPosition.setY(
                ( G4double( rowLeft ) - G4double( nCrystalsInColumn ) / 2 ) *
                                        crystalHeight + crystalHeight / 2 );
            calorimeterEPRightPosition.setX(
                ( G4double( columnRight ) - G4double( nCrystalsInRow ) / 2 ) *
                                        crystalWidth  + crystalWidth / 2 );
            calorimeterEPRightPosition.setY(
                ( G4double( rowRight ) - G4double( nCrystalsInColumn ) / 2 ) *
                                        crystalHeight + crystalHeight / 2 );
        }
        break;
    case CexmcEntryPointByLinearEDWeights :
    case CexmcEntryPointBySqrtEDWeights :
        {
            G4double  x( 0 );
            G4double  y( 0 );

            CalculateWeightedEPPosition( edStore->calorimeterEDLeftCollection,
                                         rowLeft, columnLeft, x, y,
                                         calorimeterEDLeftAdjacent );
            calorimeterEPLeftPosition.setX( x );
            calorimeterEPLeftPosition.setY( y );
            CalculateWeightedEPPosition( edStore->calorimeterEDRightCollection,
                                         rowRight, columnRight, x, y,
                                         calorimeterEDRightAdjacent );
            calorimeterEPRightPosition.setX( x );
            calorimeterEPRightPosition.setY( y );

            if ( csAlgorithm == CexmcSelectAdjacentCrystals )
                edInAdjacentCrystalsCollected = true;
        }
        break;
    default :
        break;
    }

    switch ( epDepthDefinitionAlgorithm )
    {
    case CexmcEntryPointDepthPlain :
        break;
    case CexmcEntryPointDepthSphere :
        {
            G4double  calorimeterEPLeftRadiusOfTheSphere(
                             calorimeterLeftTransform.NetTranslation().mag() +
                             calorimeterEPLeftPosition.z() );
            G4double  calorimeterEPLeftRadiusOfTheSphere2(
                                      calorimeterEPLeftRadiusOfTheSphere *
                                      calorimeterEPLeftRadiusOfTheSphere );
            G4double  calorimeterEPLeftPositionX2(
                                            calorimeterEPLeftPosition.x() *
                                            calorimeterEPLeftPosition.x() );
            G4double  calorimeterEPLeftPositionY2(
                                            calorimeterEPLeftPosition.y() *
                                            calorimeterEPLeftPosition.y() );
            G4double  calorimeterEPLeftPositionZOffset(
                           calorimeterEPLeftRadiusOfTheSphere - std::sqrt(
                                  calorimeterEPLeftRadiusOfTheSphere2 -
                                  calorimeterEPLeftPositionX2 -
                                  calorimeterEPLeftPositionY2 ) );
            G4double  calorimeterEPRightRadiusOfTheSphere(
                              calorimeterRightTransform.NetTranslation().mag() +
                              calorimeterEPRightPosition.z() );
            G4double  calorimeterEPRightRadiusOfTheSphere2(
                                       calorimeterEPRightRadiusOfTheSphere *
                                       calorimeterEPRightRadiusOfTheSphere );
            G4double  calorimeterEPRightPositionX2(
                                            calorimeterEPRightPosition.x() *
                                            calorimeterEPRightPosition.x() );
            G4double  calorimeterEPRightPositionY2(
                                            calorimeterEPRightPosition.y() *
                                            calorimeterEPRightPosition.y() );
            G4double  calorimeterEPRightPositionZOffset(
                            calorimeterEPRightRadiusOfTheSphere - std::sqrt(
                                    calorimeterEPRightRadiusOfTheSphere2 -
                                    calorimeterEPRightPositionX2 -
                                    calorimeterEPRightPositionY2 ) );
            calorimeterEPLeftPosition.setZ( calorimeterEPLeftPosition.z() -
                                            calorimeterEPLeftPositionZOffset );
            calorimeterEPRightPosition.setZ( calorimeterEPRightPosition.z() -
                                         calorimeterEPRightPositionZOffset );
        }
        break;
    default :
        break;
    }

    calorimeterEPLeftWorldPosition = calorimeterLeftTransform.TransformPoint(
                                                 calorimeterEPLeftPosition );
    calorimeterEPLeftWorldDirection = calorimeterLeftTransform.TransformAxis(
                                                 calorimeterEPLeftDirection );
    calorimeterEPRightWorldPosition = calorimeterRightTransform.TransformPoint(
                                                 calorimeterEPRightPosition );
    calorimeterEPRightWorldDirection = calorimeterRightTransform.TransformAxis(
                                                 calorimeterEPRightDirection );

    if ( collectEDInAdjacentCrystals && ! edInAdjacentCrystalsCollected )
    {
        CollectEDInAdjacentCrystals( edStore->calorimeterEDLeftCollection,
                                     rowLeft, columnLeft,
                                     calorimeterEDLeftAdjacent );
        CollectEDInAdjacentCrystals( edStore->calorimeterEDRightCollection,
                                     rowRight, columnRight,
                                     calorimeterEDRightAdjacent );
    }

    hasBasicTrigger = true;
}


void  CexmcReconstructor::ReconstructTargetPoint( void )
{
    if ( ! targetEPInitialized )
    {
        targetEPWorldPosition = targetTransform.TransformPoint(
                                                    G4ThreeVector( 0, 0, 0 ) );
        targetEPWorldDirection.setX( 0 );
        targetEPWorldDirection.setY( 0 );
        targetEPWorldDirection.setZ( 1 );

        targetEPPosition = targetTransform.Inverse().TransformPoint(
                                                    targetEPWorldPosition );
        targetEPDirection = targetTransform.Inverse().TransformAxis(
                                                    targetEPWorldDirection );
        targetEPInitialized = true;
    }

    hasBasicTrigger = true;
}


void  CexmcReconstructor::ReconstructAngle( void )
{
    G4ThreeVector  epLeft( calorimeterEPLeftWorldPosition -
                           targetEPWorldPosition );
    G4ThreeVector  epRight( calorimeterEPRightWorldPosition -
                            targetEPWorldPosition );
    theAngle = epLeft.angle( epRight );

    hasBasicTrigger = true;
}


void  CexmcReconstructor::CollectEDInAdjacentCrystals(
                        const CexmcEnergyDepositCalorimeterCollection &  edHits,
                        G4int  row, G4int  column, G4double &  ed )
{
    G4int  i( 0 );

    for ( CexmcEnergyDepositCalorimeterCollection::const_iterator
                                k( edHits.begin() ); k != edHits.end(); ++k )
    {
        if ( i - row > 1 || i - row < -1 )
        {
            ++i;
            continue;
        }

        G4int  j( 0 );
        for ( CexmcEnergyDepositCrystalRowCollection::const_iterator
                  l( k->begin() ); l != k->end(); ++l )
        {
            if ( j - column > 1 || j - column < -1 )
            {
                ++j;
                continue;
            }
            ed += *l;
            ++j;
        }
        ++i;
    }
}


void  CexmcReconstructor::CalculateWeightedEPPosition(
                const CexmcEnergyDepositCalorimeterCollection &  edHits,
                G4int  row, G4int  column, G4double &  x, G4double &  y,
                G4double &  ed )
{
    G4int     nCrystalsInColumn( calorimeterGeometry.nCrystalsInColumn );
    G4int     nCrystalsInRow( calorimeterGeometry.nCrystalsInRow );
    G4double  crystalWidth( calorimeterGeometry.crystalWidth );
    G4double  crystalHeight( calorimeterGeometry.crystalHeight );

    G4int     i( 0 );
    G4double  xWeightsSum( 0 );
    G4double  yWeightsSum( 0 );
    G4double  energyWeightsSum( 0 );

    if ( csAlgorithm == CexmcSelectAdjacentCrystals )
        ed = 0.;

    for ( CexmcEnergyDepositCalorimeterCollection::const_iterator
                                k( edHits.begin() ); k != edHits.end(); ++k )
    {
        if ( csAlgorithm == CexmcSelectAdjacentCrystals &&
                                         ( i - row > 1 || i - row < -1 ) )
        {
            ++i;
            continue;
        }

        G4int  j( 0 );
        for ( CexmcEnergyDepositCrystalRowCollection::const_iterator
                  l( k->begin() ); l != k->end(); ++l )
        {
            if ( csAlgorithm == CexmcSelectAdjacentCrystals &&
                                         ( j - column > 1 || j - column < -1 ) )
            {
                ++j;
                continue;
            }

            if ( csAlgorithm == CexmcSelectAdjacentCrystals )
                ed += *l;
            
            G4double  xInCalorimeterOffset(
                        ( G4double( j ) - G4double( nCrystalsInRow ) / 2 ) *
                        crystalWidth  + crystalWidth / 2 );
            G4double  energyWeight(
                        epDefinitionAlgorithm ==
                                        CexmcEntryPointBySqrtEDWeights ?
                                                        std::sqrt( *l ) : *l );
            xWeightsSum += energyWeight * xInCalorimeterOffset;
            G4double  yInCalorimeterOffset(
                        ( G4double( i ) - G4double( nCrystalsInColumn ) / 2 ) *
                        crystalHeight  + crystalHeight / 2 );
            yWeightsSum += energyWeight * yInCalorimeterOffset;
            energyWeightsSum += energyWeight;
            ++j;
        }
        ++i;
    }

    x = xWeightsSum / energyWeightsSum;
    y = yWeightsSum / energyWeightsSum;
}


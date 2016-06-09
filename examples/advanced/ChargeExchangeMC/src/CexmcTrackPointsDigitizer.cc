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
 *       Filename:  CexmcTrackPointsDigitizer.cc
 *
 *    Description:  track points collector
 *
 *        Version:  1.0
 *        Created:  24.11.2009 16:34:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4DigiManager.hh>
#include <G4RunManager.hh>
#include <G4String.hh>
#include "CexmcTrackPointsDigitizer.hh"
#include "CexmcTrackPoints.hh"
#include "CexmcTrackPointsInLeftRightSet.hh"
#include "CexmcTrackPointsInCalorimeter.hh"
#include "CexmcSensitiveDetectorsAttributes.hh"
#include "CexmcCommon.hh"


CexmcTrackPointsDigitizer::CexmcTrackPointsDigitizer( const G4String &  name ) :
    G4VDigitizerModule( name ), hasTriggered( false )
{
    G4RunManager *      runManager( G4RunManager::GetRunManager() );
    const CexmcSetup *  setup( static_cast< const CexmcSetup * >(
                                runManager->GetUserDetectorConstruction() ) );
    calorimeterGeometry = setup->GetCalorimeterGeometry();
}


void  CexmcTrackPointsDigitizer::InitializeData( void )
{
    monitorTP.trackId = CexmcInvalidTrackId;
    targetTPBeamParticle.trackId = CexmcInvalidTrackId;
    targetTPOutputParticle.trackId = CexmcInvalidTrackId;
    targetTPNucleusParticle.trackId = CexmcInvalidTrackId;
    targetTPOutputParticleDecayProductParticle[ 0 ].trackId =
                                                        CexmcInvalidTrackId;
    targetTPOutputParticleDecayProductParticle[ 1 ].trackId =
                                                        CexmcInvalidTrackId;
    vetoCounterTPLeft.trackId = CexmcInvalidTrackId;
    vetoCounterTPRight.trackId = CexmcInvalidTrackId;
    calorimeterTPLeft.trackId = CexmcInvalidTrackId;
    calorimeterTPRight.trackId = CexmcInvalidTrackId;
    hasTriggered = false;
}


void  CexmcTrackPointsDigitizer::Digitize( void )
{
    InitializeData();

    G4int     nCrystalsInColumn( calorimeterGeometry.nCrystalsInColumn );
    G4int     nCrystalsInRow( calorimeterGeometry.nCrystalsInRow );
    G4double  crystalWidth( calorimeterGeometry.crystalWidth );
    G4double  crystalHeight( calorimeterGeometry.crystalHeight );

    G4DigiManager *  digiManager( G4DigiManager::GetDMpointer() );
    G4int    hcId( digiManager->GetHitsCollectionID(
                    CexmcDetectorRoleName[ CexmcMonitorDetectorRole ] +
                    "/" + CexmcDetectorTypeName[ CexmcTPDetector ] ) );
    const CexmcTrackPointsCollection *
             hitsCollection( static_cast< const CexmcTrackPointsCollection * >(
                                    digiManager->GetHitsCollection( hcId ) ) );

    if ( hitsCollection )
    {
        for ( CexmcTrackPointsCollectionData::iterator
                  k( hitsCollection->GetMap()->begin() );
                      k != hitsCollection->GetMap()->end(); ++k )
        {
            monitorTP = *k->second;
            break;
        }
    }

    hcId = digiManager->GetHitsCollectionID(
                    CexmcDetectorRoleName[ CexmcTargetDetectorRole ] +
                    "/" + CexmcDetectorTypeName[ CexmcTPDetector ] );
    hitsCollection = static_cast< const CexmcTrackPointsCollection * >(
                                    digiManager->GetHitsCollection( hcId ) );

    if ( hitsCollection )
    {
        for ( CexmcTrackPointsCollectionData::iterator
                  k( hitsCollection->GetMap()->begin() );
                      k != hitsCollection->GetMap()->end(); ++k )
        {
            do
            {
                if ( k->second->trackType == CexmcBeamParticleTrack )
                {
                    targetTPBeamParticle = *k->second;
                    break;
                }
                if ( k->second->trackType == CexmcOutputParticleTrack )
                {
                    targetTPOutputParticle = *k->second;
                    hasTriggered = targetTPOutputParticle.IsValid();
                    break;
                }
                if ( k->second->trackType == CexmcNucleusParticleTrack )
                {
                    targetTPNucleusParticle = *k->second;
                    break;
                }
                if ( k->second->trackType ==
                                        CexmcOutputParticleDecayProductTrack )
                {
                    /* NB: if there are more than 2 output particle's decay
                     * products then the chosen particles may differ from those
                     * which entered calorimeters; however this is not a
                     * critical issue as far as information about decay products
                     * is not necessary in reconstruction and only used in some
                     * histograming */
                    G4int  index(
                            targetTPOutputParticleDecayProductParticle[ 0 ].
                                                          trackId > 0 ? 1 : 0 );
                    targetTPOutputParticleDecayProductParticle[ index ] =
                            *k->second;
                    break;
                }
            } while ( false );
        }
    }

    hcId = digiManager->GetHitsCollectionID(
                    CexmcDetectorRoleName[ CexmcVetoCounterDetectorRole ] +
                    "/" + CexmcDetectorTypeName[ CexmcTPDetector ] );
    hitsCollection = static_cast< const CexmcTrackPointsCollection * >(
                                    digiManager->GetHitsCollection( hcId ) );

    if ( hitsCollection )
    {
        for ( CexmcTrackPointsCollectionData::iterator
                  k( hitsCollection->GetMap()->begin() );
                      k != hitsCollection->GetMap()->end(); ++k )
        {
            if ( k->second->trackType != CexmcOutputParticleDecayProductTrack )
                continue;

            G4int  index( k->first );
            CexmcSide  side( CexmcTrackPointsInLeftRightSet::GetSide(
                                                                   index ) );
            switch ( side )
            {
            case CexmcLeft :
                vetoCounterTPLeft = *k->second;
                break;
            case CexmcRight :
                vetoCounterTPRight = *k->second;
                break;
            default :
                break;
            }
        }
    }

    hcId = digiManager->GetHitsCollectionID(
                    CexmcDetectorRoleName[ CexmcCalorimeterDetectorRole ] +
                    "/" + CexmcDetectorTypeName[ CexmcTPDetector ] );
    hitsCollection = static_cast< const CexmcTrackPointsCollection * >(
                                    digiManager->GetHitsCollection( hcId ) );

    if ( hitsCollection )
    {
        for ( CexmcTrackPointsCollectionData::iterator
                  k( hitsCollection->GetMap()->begin() );
                      k != hitsCollection->GetMap()->end(); ++k )
        {
            if ( k->second->trackType != CexmcOutputParticleDecayProductTrack )
                continue;

            G4int  index( k->first );
            CexmcSide  side( CexmcTrackPointsInLeftRightSet::GetSide(
                                                                   index ) );
            G4int      row( CexmcTrackPointsInCalorimeter::GetRow( index ) );
            G4int      column( CexmcTrackPointsInCalorimeter::GetColumn(
                                                                   index ) );
            G4double   xInCalorimeterOffset(
                    ( G4double( column ) - G4double( nCrystalsInRow ) / 2 ) *
                                        crystalWidth + crystalWidth / 2 );
            G4double   yInCalorimeterOffset(
                    ( G4double( row ) - G4double( nCrystalsInColumn ) / 2 ) *
                                        crystalHeight + crystalHeight / 2 );
            switch ( side )
            {
            case CexmcLeft :
                calorimeterTPLeft = *k->second;
                calorimeterTPLeft.positionLocal.setX( xInCalorimeterOffset +
                                        calorimeterTPLeft.positionLocal.x() );
                calorimeterTPLeft.positionLocal.setY( yInCalorimeterOffset +
                                        calorimeterTPLeft.positionLocal.y() );
                break;
            case CexmcRight :
                calorimeterTPRight = *k->second;
                calorimeterTPRight.positionLocal.setX( xInCalorimeterOffset +
                                        calorimeterTPRight.positionLocal.x() );
                calorimeterTPRight.positionLocal.setY( yInCalorimeterOffset +
                                        calorimeterTPRight.positionLocal.y() );
                break;
            default :
                break;
            }
        }
    }
}


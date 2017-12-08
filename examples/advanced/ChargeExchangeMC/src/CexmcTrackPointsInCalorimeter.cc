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
 *       Filename:  CexmcTrackPointsInCalorimeter.cc
 *
 *    Description:  track points in calorimeter
 *
 *        Version:  1.0
 *        Created:  22.11.2009 21:59:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4VTouchable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4NavigationHistory.hh>
#include <G4UnitsTable.hh>
#include "CexmcTrackPointsInCalorimeter.hh"
#include "CexmcSetup.hh"


CexmcTrackPointsInCalorimeter::CexmcTrackPointsInCalorimeter(
                        const G4String &  name, const CexmcSetup *  setup_ ) :
    CexmcTrackPointsInLeftRightSet( name, setup_ )
{
}


G4int  CexmcTrackPointsInCalorimeter::GetIndex( G4Step *  step )
{
    G4int                        ret( GetTrackId( step ) );
    G4StepPoint *                preStep( step->GetPreStepPoint() );
    const G4VTouchable *         touchable( preStep->GetTouchable() );
    const G4NavigationHistory *  navHistory( touchable->GetHistory() );
    G4int                        navDepth( navHistory->GetDepth() );
    G4VPhysicalVolume *          pVolume( navHistory->GetVolume(
                                                        navDepth - 2 ) );

    if ( setup->IsRightCalorimeter( pVolume ) )
        ret |= 1 << leftRightBitsOffset;

    ret |= touchable->GetReplicaNumber( 0 ) << copyDepth0BitsOffset;
    ret |= touchable->GetReplicaNumber( 1 ) << copyDepth1BitsOffset;

    return ret;
}


void  CexmcTrackPointsInCalorimeter::PrintAll( void )
{
    G4int   nmbOfEntries( eventMap->entries() );

    if ( nmbOfEntries == 0 )
        return;

    PrintHeader( nmbOfEntries );

    for ( CexmcTrackPointsCollectionData::iterator
                         itr( eventMap->GetMap()->begin() );
                                     itr != eventMap->GetMap()->end(); ++itr )
    {
        G4bool  isRightDetector( itr->first >> leftRightBitsOffset );
        G4int   index( itr->first &
                            ( ( 1 << ( leftRightBitsOffset - 1 ) ) |
                              ( ( 1 << ( leftRightBitsOffset - 1 ) ) - 1 ) ) );
        G4int   copyDepth1( index >> copyDepth1BitsOffset );
        index &= ( 1 << ( copyDepth1BitsOffset - 1 ) ) |
                   ( ( 1 << ( copyDepth1BitsOffset - 1 ) ) - 1 );
        G4int   copyDepth0( index >> copyDepth0BitsOffset );
        G4int   trackId( index &
                         ( ( 1 << ( copyDepth0BitsOffset - 1 ) ) |
                           ( ( 1 << ( copyDepth0BitsOffset - 1 ) ) - 1 ) ) );
        const G4String  detectorSide( isRightDetector ? "right" : "left" );
        G4cout << "       " << detectorSide << " detector, row " <<
                copyDepth1 << ", column " << copyDepth0 << G4endl;
        G4cout << "         , track id " << trackId << G4endl;
        G4cout << "         , position: " <<
                G4BestUnit( itr->second->positionLocal, "Length" ) << G4endl;
        G4cout << "         , direction: " <<
                itr->second->directionLocal << G4endl;
        G4cout << "         , momentum: " <<
                G4BestUnit( itr->second->momentumAmp, "Energy" ) << G4endl;
        G4cout << "         , particle: "
                << itr->second->particle->GetParticleName() << G4endl;
    }
}


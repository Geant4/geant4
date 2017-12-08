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
 *       Filename:  CexmcRunAction.cc
 *
 *    Description:  run action
 *
 *        Version:  1.0
 *        Created:  20.12.2009 00:18:05
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <limits>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "CexmcRunAction.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcProductionModel.hh"
#include "CexmcAngularRange.hh"
#include "CexmcException.hh"


CexmcRunAction::CexmcRunAction( CexmcPhysicsManager *  physicsManager_ ) :
    physicsManager( physicsManager_ )
{
}


G4Run *  CexmcRunAction::GenerateRun( void )
{
    return new CexmcRun;
}


void  CexmcRunAction::PrintResults(
                    const CexmcNmbOfHitsInRanges &  nmbOfHitsSampled,
                    const CexmcNmbOfHitsInRanges &  nmbOfHitsSampledFull,
                    const CexmcNmbOfHitsInRanges &  nmbOfHitsTriggeredRealRange,
                    const CexmcNmbOfHitsInRanges &  nmbOfHitsTriggeredRecRange,
                    const CexmcNmbOfHitsInRanges &  nmbOfOrphanHits,
                    const CexmcAngularRangeList &  angularRanges,
                    G4int  nmbOfFalseHitsTriggeredEDT,
                    G4int  nmbOfFalseHitsTriggeredRec )
{
    /* there are 7 auxiliary columns:
     * 1. acc real, [floating point number]
     * 2. triggered real,
     * 3. total hits sampled and monitored,
     * 4. acc reconstructed, [floating point number]
     * 5. triggered reconstructed,
     * 6. total hits sampled and monitored (identical with #3),
     * 7. total hits sampled.
     * As far as #3 and #6 are identical, nmbOfAuxColumns = 6 */
    const size_t                               nmbOfAuxColumns( 6 );
    const std::streamsize                      prec( 8 );
    std::vector< std::vector< std::string > >  auxStrings;
    size_t                                     maxSize[ nmbOfAuxColumns ];

    for ( size_t  i( 0 ); i < nmbOfAuxColumns; ++i )
        maxSize[ i ] = 0;

    /* addition of 2 (for '0.') for acceptances, that are floating point
     * numbers, is correct as far as ios::fixed will be used, and no negative
     * values are expected, and values will be less than 1. */
    maxSize[ 0 ] = prec + 2;
    maxSize[ 3 ] = prec + 2;

    for ( CexmcAngularRangeList::const_iterator  k( angularRanges.begin() );
                                                k != angularRanges.end(); ++k )
    {
        G4int     total( 0 );
        G4int     totalFull( 0 );
        G4int     triggered( 0 );
        G4double  acc( std::numeric_limits< G4double >::quiet_NaN() );

        CexmcNmbOfHitsInRanges::const_iterator  found(
                                        nmbOfHitsSampled.find( k->index ) );
        if ( found != nmbOfHitsSampled.end() )
        {
            total = found->second;
            acc = 0;
        }

        found = nmbOfHitsSampledFull.find( k->index );
        if ( found != nmbOfHitsSampledFull.end() )
        {
            totalFull = found->second;
        }

        G4double  accSave( acc );
        found = nmbOfHitsTriggeredRealRange.find( k->index );
        if ( found != nmbOfHitsTriggeredRealRange.end() )
        {
            triggered = found->second;
            if ( total > 0 )
                acc = G4double( triggered ) / total;
        }

        std::ostringstream  auxStringStream[ nmbOfAuxColumns ];

        for ( size_t  i( 0 ); i < nmbOfAuxColumns; ++i )
        {
            auxStringStream[ i ].precision( prec );
            auxStringStream[ i ].flags( std::ios::fixed );
        }

        G4int  i( 0 );

        auxStringStream[ i ] << acc;
        auxStringStream[ ++i ] << triggered;
        size_t  size( auxStringStream[ i ].str().size() );
        maxSize[ i ] = maxSize[ i ] > size ? maxSize[ i ] : size;
        auxStringStream[ ++i ] << total;
        size = auxStringStream[ i ].str().size();
        maxSize[ i ] = maxSize[ i ] > size ? maxSize[ i ] : size;

        triggered = 0;
        acc = accSave;
        found = nmbOfHitsTriggeredRecRange.find( k->index );
        if ( found != nmbOfHitsTriggeredRecRange.end() )
        {
            triggered = found->second;
            if ( total > 0 )
                acc = G4double( triggered ) / total;
        }

        auxStringStream[ ++i ] << acc;
        auxStringStream[ ++i ] << triggered;
        size = auxStringStream[ i ].str().size();
        maxSize[ i ] = maxSize[ i ] > size ? maxSize[ i ] : size;
        auxStringStream[ ++i ] << totalFull;
        size = auxStringStream[ i ].str().size();
        maxSize[ i ] = maxSize[ i ] > size ? maxSize[ i ] : size;

        std::vector< std::string >  auxString( nmbOfAuxColumns );

        for ( size_t  j( 0 ); j < nmbOfAuxColumns; ++j )
            auxString[ j ] = auxStringStream[ j ].str();

        auxStrings.push_back( auxString );
    }

    G4cout << " --- Setup acceptances (range | real (trg / mon) | "
              "rec (trg / mon / all)):" << G4endl;

    G4int  i( 0 );
    for ( CexmcAngularRangeList::const_iterator  k( angularRanges.begin() );
                                                k != angularRanges.end(); ++k )
    {
        G4cout << "       " << *k;
        G4int  j( 0 );
        G4cout << "  | " << std::setw( maxSize[ j ] );
        G4cout << auxStrings[ i ][ j++ ];
        G4cout << " ( " << std::setw( maxSize[ j ] );
        G4cout << auxStrings[ i ][ j++ ];
        G4cout << " / " << std::setw( maxSize[ j ] );
        G4cout << auxStrings[ i ][ j++ ];
        G4cout << " )  | " << std::setw( maxSize[ j ] );
        G4cout << auxStrings[ i ][ j++ ];
        G4cout << " ( " << std::setw( maxSize[ j ] );
        G4cout << auxStrings[ i ][ j++ ];
        G4cout << " / " << std::setw( maxSize[ 2 ] );
        G4cout << auxStrings[ i ][ 2 ];
        G4cout << " / " << std::setw( maxSize[ j ] );
        G4cout << auxStrings[ i++ ][ j++ ] << " )" << G4endl;
    }

    CexmcAngularRangeList  angularGaps;
    GetAngularGaps( angularRanges, angularGaps );

    if ( ! angularGaps.empty() )
    {
        G4cout << "    orphans detected: " << G4endl;
        for ( CexmcAngularRangeList::const_iterator  k( angularGaps.begin() );
                                                k != angularGaps.end(); ++k )
        {
            G4cout << "       " << *k;
            G4int     total( 0 );

            CexmcNmbOfHitsInRanges::const_iterator  found(
                                            nmbOfOrphanHits.find( k->index ) );
            if ( found != nmbOfHitsSampled.end() )
            {
                total = found->second;
            }
            G4cout << " " << total << G4endl;
        }
    }

    G4cout << "       ---" << G4endl;
    G4cout << "       False hits (edt | rec):  " <<
        nmbOfFalseHitsTriggeredEDT << " | " << nmbOfFalseHitsTriggeredRec <<
        G4endl;
}


void  CexmcRunAction::EndOfRunAction( const G4Run *  run )
{
    const CexmcRun *  theRun( static_cast< const CexmcRun * >( run ) );

    const CexmcNmbOfHitsInRanges &  nmbOfHitsSampled(
                                        theRun->GetNmbOfHitsSampled() );
    const CexmcNmbOfHitsInRanges &  nmbOfHitsSampledFull(
                                        theRun->GetNmbOfHitsSampledFull() );
    const CexmcNmbOfHitsInRanges &  nmbOfHitsTriggeredRealRange(
                                    theRun->GetNmbOfHitsTriggeredRealRange() );
    const CexmcNmbOfHitsInRanges &  nmbOfHitsTriggeredRecRange(
                                    theRun->GetNmbOfHitsTriggeredRecRange() );
    const CexmcNmbOfHitsInRanges &  nmbOfOrphanHits(
                                        theRun->GetNmbOfOrphanHits() );

    CexmcProductionModel *          productionModel(
                                        physicsManager->GetProductionModel() );
    if ( ! productionModel )
            throw CexmcException( CexmcWeirdException );

    const CexmcAngularRangeList &  angularRanges(
                                        productionModel->GetAngularRanges() );

    G4cout << G4endl;
    PrintResults( nmbOfHitsSampled, nmbOfHitsSampledFull,
                  nmbOfHitsTriggeredRealRange, nmbOfHitsTriggeredRecRange,
                  nmbOfOrphanHits, angularRanges,
                  theRun->GetNmbOfFalseHitsTriggeredEDT(),
                  theRun->GetNmbOfFalseHitsTriggeredRec() );
    G4cout << G4endl;
}


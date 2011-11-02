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
 *       Filename:  CexmcRun.cc
 *
 *    Description:  run data (acceptances etc.)
 *
 *        Version:  1.0
 *        Created:  19.12.2009 23:59:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcRun.hh"


CexmcRun::CexmcRun() : nmbOfFalseHitsTriggeredEDT( 0 ),
    nmbOfFalseHitsTriggeredRec( 0 ), nmbOfSavedEvents( 0 ),
    nmbOfSavedFastEvents( 0 )
{
}


void  CexmcRun::IncrementNmbOfHitsSampled( G4int  index )
{
    CexmcNmbOfHitsInRanges::iterator  found(
                                        nmbOfHitsSampled.find( index ) );
    if ( found == nmbOfHitsSampled.end() )
        nmbOfHitsSampled.insert( CexmcNmbOfHitsInRangesData( index, 1 ) );
    else
        ++found->second;
}


void  CexmcRun::IncrementNmbOfHitsSampledFull( G4int  index )
{
    CexmcNmbOfHitsInRanges::iterator  found(
                                        nmbOfHitsSampledFull.find( index ) );
    if ( found == nmbOfHitsSampledFull.end() )
        nmbOfHitsSampledFull.insert( CexmcNmbOfHitsInRangesData( index, 1 ) );
    else
        ++found->second;
}


void  CexmcRun::IncrementNmbOfHitsTriggeredRealRange( G4int  index )
{
    CexmcNmbOfHitsInRanges::iterator  found(
                                    nmbOfHitsTriggeredRealRange.find( index ) );
    if ( found == nmbOfHitsTriggeredRealRange.end() )
        nmbOfHitsTriggeredRealRange.insert(
                                    CexmcNmbOfHitsInRangesData( index, 1 ) );
    else
        ++found->second;
}


void  CexmcRun::IncrementNmbOfHitsTriggeredRecRange( G4int  index )
{
    CexmcNmbOfHitsInRanges::iterator  found(
                                    nmbOfHitsTriggeredRecRange.find( index ) );
    if ( found == nmbOfHitsTriggeredRecRange.end() )
        nmbOfHitsTriggeredRecRange.insert(
                                    CexmcNmbOfHitsInRangesData( index, 1 ) );
    else
        ++found->second;
}


void  CexmcRun::IncrementNmbOfOrphanHits( G4int  index )
{
    CexmcNmbOfHitsInRanges::iterator  found(
                                        nmbOfOrphanHits.find( index ) );
    if ( found == nmbOfOrphanHits.end() )
        nmbOfOrphanHits.insert( CexmcNmbOfHitsInRangesData( index, 1 ) );
    else
        ++found->second;
}


void  CexmcRun::IncrementNmbOfFalseHitsTriggeredEDT( void )
{
    ++nmbOfFalseHitsTriggeredEDT;
}


void  CexmcRun::IncrementNmbOfFalseHitsTriggeredRec( void )
{
    ++nmbOfFalseHitsTriggeredRec;
}


void  CexmcRun::IncrementNmbOfSavedEvents( void )
{
    ++nmbOfSavedEvents;
}


void  CexmcRun::IncrementNmbOfSavedFastEvents( void )
{
    ++nmbOfSavedFastEvents;
}


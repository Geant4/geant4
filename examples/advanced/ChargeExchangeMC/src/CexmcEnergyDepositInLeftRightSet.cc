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
 *       Filename:  CexmcEnergyDepositInLeftRightSet.cc
 *
 *    Description:  energy deposit scorer in left/right detector sets
 *                  (e.g. veto counters and calorimeters)
 *
 *        Version:  1.0
 *        Created:  14.11.2009 12:48:22
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
#include <G4UnitsTable.hh>
#include "CexmcEnergyDepositInLeftRightSet.hh"
#include "CexmcSetup.hh"


CexmcEnergyDepositInLeftRightSet::CexmcEnergyDepositInLeftRightSet(
                        const G4String &  name, const CexmcSetup *  setup_ ) :
    CexmcSimpleEnergyDeposit( name ), setup( setup_ )
{
}


G4int  CexmcEnergyDepositInLeftRightSet::GetIndex( G4Step *  step )
{
    G4int                ret( 0 );
    G4StepPoint *        preStep( step->GetPreStepPoint() );
    G4VPhysicalVolume *  pVolume( preStep->GetPhysicalVolume() );

    if ( setup->IsRightDetector( pVolume ) )
        ret |= 1 << leftRightBitsOffset;

    return ret;
}


void  CexmcEnergyDepositInLeftRightSet::PrintAll( void )
{
    G4int   nmbOfEntries( eventMap->entries() );

    if ( nmbOfEntries == 0 )
        return;

    PrintHeader( nmbOfEntries );

    for ( CexmcEnergyDepositCollectionData::iterator
                         itr( eventMap->GetMap()->begin() );
                                     itr != eventMap->GetMap()->end(); ++itr )
    {
        G4bool  isRightDetector( itr->first >> leftRightBitsOffset );
        const G4String  detectorSide( isRightDetector ? "right" : "left" );
        G4cout << "       " << detectorSide << " detector" << G4endl;
        G4cout << "         , energy deposit " <<
                G4BestUnit( *( itr->second ), "Energy" ) << G4endl;
    }
}


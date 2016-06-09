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
 *       Filename:  CexmcSimpleEnergyDeposit.cc
 *
 *    Description:  simple energy deposit scorer
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

#include <G4UnitsTable.hh>
#include <G4HCofThisEvent.hh>
#include "CexmcSimpleEnergyDeposit.hh"


CexmcSimpleEnergyDeposit::CexmcSimpleEnergyDeposit( const G4String &  name ) :
    CexmcPrimitiveScorer( name ), hcId( -1 )
{
}


G4int  CexmcSimpleEnergyDeposit::GetIndex( G4Step * )
{
    return 0;
}


G4bool  CexmcSimpleEnergyDeposit::ProcessHits( G4Step *  step,
                                               G4TouchableHistory * )
{
    G4double  energyDeposit( step->GetTotalEnergyDeposit() );

    if ( energyDeposit == 0. )
        return false;

    eventMap->add( GetIndex( step ), energyDeposit ); 

    return true; 
}


void  CexmcSimpleEnergyDeposit::Initialize( G4HCofThisEvent *  hcOfEvent )
{
    eventMap = new CexmcEnergyDepositCollection( detector->GetName(),
                                                 primitiveName );

    if ( hcId < 0 )
        hcId = GetCollectionID( 0 );

    hcOfEvent->AddHitsCollection( hcId, eventMap );
}


void  CexmcSimpleEnergyDeposit::EndOfEvent( G4HCofThisEvent * )
{
    if ( GetVerboseLevel() > 0 )
        PrintAll();
}


void  CexmcSimpleEnergyDeposit::clear( void )
{
    eventMap->clear();
}


void  CexmcSimpleEnergyDeposit::DrawAll( void )
{
}


void  CexmcSimpleEnergyDeposit::PrintAll( void )
{
    G4int   nmbOfEntries( eventMap->entries() );

    if ( nmbOfEntries == 0 )
        return;

    PrintHeader( nmbOfEntries );

    /* index is always 0 */
    for ( CexmcEnergyDepositCollectionData::iterator
                         itr( eventMap->GetMap()->begin() );
                                     itr != eventMap->GetMap()->end(); ++itr )
    {
        G4cout << "       energy deposit: " <<
                G4BestUnit( *( itr->second ), "Energy" ) << G4endl;
    }
}


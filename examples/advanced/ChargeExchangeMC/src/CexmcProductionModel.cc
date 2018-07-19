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
 *       Filename:  CexmcProductionModel.cc
 *
 *    Description:  interface to production model
 *
 *        Version:  1.0
 *        Created:  03.11.2009 16:56:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcRunManager.hh"
#include "CexmcProductionModel.hh"
#include "CexmcProductionModelMessenger.hh"


CexmcProductionModel::CexmcProductionModel( const G4String &  name_,
                                            G4bool  fermiMotionIsOn_ ) :
    name( name_ ), fermiMotionIsOn( fermiMotionIsOn_ ),
    incidentParticle( NULL ), nucleusParticle( NULL ), outputParticle( NULL ),
    nucleusOutputParticle( NULL ), messenger( NULL )
{
    angularRanges.push_back( CexmcAngularRange( 1.0, -1.0, 0 ) );
    messenger = new CexmcProductionModelMessenger( this );
}


CexmcProductionModel::~CexmcProductionModel()
{
    delete messenger;
}


void  CexmcProductionModel::SetAngularRange( G4double  top, G4double  bottom,
                                             G4int  nmbOfDivs )
{
    if ( ! IsValidCandidateForAngularRange( top, bottom, nmbOfDivs ) )
        throw CexmcException( CexmcInvalidAngularRange );

    if ( ! IsGoodCandidateForAngularRange( top, bottom ) )
        throw CexmcException( CexmcBadAngularRange );

    angularRanges.clear();
    G4double  curBottom( top );
    for ( int  i( 0 ); i < nmbOfDivs; ++i )
    {
        G4double  binWidth( ( top - bottom ) / nmbOfDivs );
        G4double  curTop( curBottom );
        curBottom -= binWidth;
        angularRanges.push_back( CexmcAngularRange( curTop, curBottom, i ) );
    }
#ifdef CEXMC_USE_ROOT
    CexmcHistoManager::Instance()->SetupARHistos( angularRanges );
#endif
}


void  CexmcProductionModel::AddAngularRange( G4double  top, G4double  bottom,
                                             G4int  nmbOfDivs )
{
    if ( ! IsValidCandidateForAngularRange( top, bottom, nmbOfDivs ) )
        throw CexmcException( CexmcInvalidAngularRange );

    if ( ! IsGoodCandidateForAngularRange( top, bottom ) )
        throw CexmcException( CexmcBadAngularRange );

    G4int  curIndex( angularRanges.size() );
    G4double  curBottom( top );
    for ( int  i( 0 ); i < nmbOfDivs; ++i )
    {
        G4double  binWidth( ( top - bottom ) / nmbOfDivs );
        G4double  curTop( curBottom );
        curBottom -= binWidth;
        CexmcAngularRange  aRange( curTop, curBottom, curIndex + i );
        angularRanges.push_back( aRange );
#ifdef CEXMC_USE_ROOT
        CexmcHistoManager::Instance()->AddARHistos( aRange );
#endif
    }
}


void  CexmcProductionModel::SetTriggeredAngularRanges( G4double  opCosThetaSCM )
{
    triggeredAngularRanges.clear();

    for ( CexmcAngularRangeList::iterator  k( angularRanges.begin() );
                                            k != angularRanges.end(); ++k )
    {
        if ( opCosThetaSCM <= k->top && opCosThetaSCM > k->bottom )
            triggeredAngularRanges.push_back( CexmcAngularRange(
                                            k->top, k->bottom, k->index ) );
    }
}


void  CexmcProductionModel::FermiMotionStatusChangeHook( void )
{
}


G4bool  CexmcProductionModel::IsGoodCandidateForAngularRange( G4double  top,
                                                      G4double  bottom ) const
{
    CexmcRunManager *  runManager( static_cast< CexmcRunManager * >(
                                           G4RunManager::GetRunManager() ) );

    if ( ! runManager->ProjectIsRead() )
        return true;

    CexmcAngularRangeList  normalizedARanges;
    GetNormalizedAngularRange( angularRangesRef, normalizedARanges );

    for ( CexmcAngularRangeList::iterator  k( normalizedARanges.begin() );
                                            k != normalizedARanges.end(); ++k )
    {
        if ( top <= k->top && bottom >= k->bottom )
            return true;
    }

    return false;
}


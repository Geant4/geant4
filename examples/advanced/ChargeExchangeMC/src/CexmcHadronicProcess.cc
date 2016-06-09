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
 *       Filename:  CexmcHadronicProcess.cc
 *
 *    Description:  hadronic process with production model
 *
 *        Version:  1.0
 *        Created:  31.10.2009 23:54:38
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4ParticleChange.hh>
#include <G4ParticleDefinition.hh>
#include <G4HadronicInteraction.hh>
#include <G4Track.hh>
#include <G4Step.hh>
#include <G4Element.hh>
#include <G4StableIsotopes.hh>
#include <G4TrackStatus.hh>
#include "CexmcHadronicProcess.hh"
#include "CexmcProductionModel.hh"
#include "CexmcIncidentParticleTrackInfo.hh"
#include "CexmcException.hh"


CexmcHadronicProcess::CexmcHadronicProcess( const G4String &  name ) :
    G4HadronicProcess( name ), productionModel( NULL ), interaction( NULL ),
    theTotalResult( NULL ), isInitialized( false )
{
    theTotalResult = new G4ParticleChange();
}


CexmcHadronicProcess::~CexmcHadronicProcess()
{
    delete theTotalResult;
}


void  CexmcHadronicProcess::RegisterProductionModel(
                                                CexmcProductionModel *  model )
{
    productionModel = model;

    interaction = dynamic_cast< G4HadronicInteraction * >( productionModel );

    if ( ! interaction )
        throw CexmcException( CexmcIncompatibleProductionModel );

    G4HadronicProcess::RegisterMe( interaction );
}


void  CexmcHadronicProcess::CalculateTargetNucleus(
                                                const G4Material *  material )
{
    G4int  numberOfElements( material->GetNumberOfElements() );
    if ( numberOfElements > 1 )
    {
        G4cout << CEXMC_LINE_START "WARNING: Number of elements in target "
                  "material is more than 1.\n              Only the first "
                  "element will be chosen for target nucleus" << G4endl;
    }

    const G4Element *  element( material->GetElement( 0 ) );
    G4double           ZZ( element->GetZ() );
    G4int              Z( G4int( ZZ + 0.5 ) );

    G4StableIsotopes  stableIsotopes;
    G4int             index( stableIsotopes.GetFirstIsotope( Z ) );
    G4double          AA( stableIsotopes.GetIsotopeNucleonCount( index ) );

    targetNucleus.SetParameters( AA, ZZ );
}


void  CexmcHadronicProcess::FillTotalResult( G4HadFinalState *  hadFinalState,
                                             const G4Track &  track )
{
    G4int  numberOfSecondaries( hadFinalState->GetNumberOfSecondaries() );

    theTotalResult->Clear();
    theTotalResult->Initialize( track );
    theTotalResult->SetSecondaryWeightByProcess( true );
    theTotalResult->ProposeLocalEnergyDeposit(
                                    hadFinalState->GetLocalEnergyDeposit() );
    theTotalResult->SetNumberOfSecondaries( numberOfSecondaries );
    theTotalResult->ProposeEnergy( hadFinalState->GetEnergyChange() );
    theTotalResult->ProposeTrackStatus( fAlive );
    if ( hadFinalState->GetStatusChange() == stopAndKill )
        theTotalResult->ProposeTrackStatus( fStopAndKill );

    for ( G4int  i( 0 ); i < numberOfSecondaries; ++i )
    {
        G4double   time( hadFinalState->GetSecondary( i )->GetTime() );
        if ( time < 0 )
            time = track.GetGlobalTime();

        G4Track *  newTrack( new G4Track(
                             hadFinalState->GetSecondary( i )->GetParticle(),
                             time, track.GetPosition() ) );

        G4double   newWeight( track.GetWeight() *
                              hadFinalState->GetSecondary( i )->GetWeight() );
        newTrack->SetWeight( newWeight );
        newTrack->SetTouchableHandle( track.GetTouchableHandle() );
        theTotalResult->AddSecondary( newTrack );
    }

    hadFinalState->Clear();
}


G4VParticleChange *  CexmcHadronicProcess::PostStepDoIt( const G4Track &  track,
                                                         const G4Step & )
{
    G4TrackStatus  trackStatus( track.GetTrackStatus() );

    if ( trackStatus != fAlive && trackStatus != fSuspend )
    {
        theTotalResult->Clear();
        theTotalResult->Initialize( track );

        return theTotalResult;
    }

    /* NB: the target nucleus is chosen only once, it means that it will always
     * have same Z and A, practically the first stable isotope of the first
     * element in elements vector will be chosen. This simplification prompts
     * the user to choose simple single-element material for the target, for
     * example liquid hydrogen. On the other hand target nucleus is supposedly
     * only needed if user decides to turn Fermi motion on, so this
     * simplification should not be very harmful */
    if ( ! isInitialized )
    {
        CalculateTargetNucleus( track.GetMaterial() );
        isInitialized = true;
    }

    G4HadProjectile    projectile( track );
    G4HadFinalState *  result( interaction->ApplyYourself( projectile,
                                                           targetNucleus ) );
    FillTotalResult( result, track );

    if ( theTotalResult->GetTrackStatus() != fStopAndKill )
    {
        CexmcTrackInfo *  trackInfo( static_cast< CexmcTrackInfo * >(
                                                track.GetUserInformation() ) );

        if ( trackInfo &&
             trackInfo->GetTypeInfo() == CexmcIncidentParticleTrackType )
        {
            CexmcIncidentParticleTrackInfo *  theTrackInfo(
                static_cast< CexmcIncidentParticleTrackInfo * >( trackInfo ) );
            theTrackInfo->SetNeedsTrackLengthResampling();
        }
    }

    return theTotalResult;
}


G4bool  CexmcHadronicProcess::IsApplicable(
                                        const G4ParticleDefinition &  particle )
{
    if ( ! productionModel )
        return false;

    G4ParticleDefinition *  incidentParticle(
                                    productionModel->GetIncidentParticle() );

    if ( ! incidentParticle )
        return false;

    return particle == *incidentParticle;
}


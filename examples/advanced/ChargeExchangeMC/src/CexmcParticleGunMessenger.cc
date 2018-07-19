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
 *       Filename:  CexmcParticleGunMessenger.cc
 *
 *    Description:  original position and momentum of the incident beam particle
 *
 *        Version:  1.0
 *        Created:  15.12.2009 14:02:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include "CexmcParticleGun.hh"
#include "CexmcParticleGunMessenger.hh"
#include "CexmcRunManager.hh"
#include "CexmcException.hh"
#include "CexmcMessenger.hh"


CexmcParticleGunMessenger::CexmcParticleGunMessenger(
                                            CexmcParticleGun *  particleGun_ ) :
    particleGun( particleGun_ ), setParticle( NULL ), setOrigPosition( NULL ),
    setOrigDirection( NULL ), setOrigMomentumAmp( NULL )
{
    setParticle = new G4UIcmdWithAString(
        ( CexmcMessenger::gunDirName + "particle" ).c_str(), this );
    setParticle->SetGuidance( "Incident beam particle" );
    setParticle->SetParameterName( "BeamParticle", false );
    setParticle->SetCandidates( "pi-" );
    setParticle->SetDefaultValue( "pi-" );
    setParticle->AvailableForStates( G4State_PreInit, G4State_Idle );

    setOrigPosition = new G4UIcmdWith3VectorAndUnit( 
        ( CexmcMessenger::gunDirName + "position" ).c_str(), this );
    setOrigPosition->SetGuidance( "Original position of the beam" );
    setOrigPosition->SetParameterName( "PositionX", "PositionY", "PositionZ",
                                       false );
    setOrigPosition->SetUnitCandidates( "mm cm m" );
    setOrigPosition->SetDefaultUnit( "cm" );
    setOrigPosition->AvailableForStates( G4State_PreInit, G4State_Idle );

    setOrigDirection = new G4UIcmdWith3Vector(
        ( CexmcMessenger::gunDirName + "direction" ).c_str(), this );
    setOrigDirection->SetGuidance( "Original direction of the beam" );
    setOrigDirection->SetParameterName( "DirectionX", "DirectionY",
                                        "DirectionZ", false );
    setOrigDirection->SetRange(
        "DirectionX >= -1.0 && DirectionX <= 1.0 && "
        "DirectionY >= -1.0 && DirectionY <= 1.0 && "
        "DirectionZ >= -1.0 && DirectionZ <= 1.0" );
    setOrigDirection->AvailableForStates( G4State_PreInit, G4State_Idle );

    setOrigMomentumAmp = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::gunDirName + "momentumAmp" ).c_str(), this );
    setOrigMomentumAmp->SetGuidance( "Original momentum of the beam" );
    setOrigMomentumAmp->SetParameterName( "MomentumAmp", false );
    setOrigMomentumAmp->SetRange( "MomentumAmp > 0" );
    setOrigMomentumAmp->SetUnitCandidates( "eV keV MeV GeV" );
    setOrigMomentumAmp->SetDefaultUnit( "MeV" );
    setOrigMomentumAmp->AvailableForStates( G4State_PreInit, G4State_Idle );
}


CexmcParticleGunMessenger::~CexmcParticleGunMessenger()
{
    delete setParticle;
    delete setOrigPosition;
    delete setOrigDirection;
    delete setOrigMomentumAmp;
}


void  CexmcParticleGunMessenger::SetNewValue( G4UIcommand *  cmd,
                                              G4String  value )
{
    do
    {
        if ( cmd == setParticle )
        {
            G4ParticleDefinition *  particleDefinition(
                G4ParticleTable::GetParticleTable()->FindParticle( value ) );

            if ( ! particleDefinition )
                throw CexmcException( CexmcWeirdException );

            particleGun->SetBeamParticle( particleDefinition );

            CexmcRunManager *  runManager( static_cast< CexmcRunManager * >(
                                            G4RunManager::GetRunManager() ) );
            runManager->BeamParticleChangeHook();
            break;
        }
        if ( cmd == setOrigPosition )
        {
            particleGun->SetOrigPosition(
                    G4UIcmdWith3VectorAndUnit::GetNew3VectorValue( value ) );
            break;
        }
        if ( cmd == setOrigDirection )
        {
            particleGun->SetOrigDirection(
                    G4UIcmdWith3Vector::GetNew3VectorValue( value ) );
            break;
        }
        if ( cmd == setOrigMomentumAmp )
        {
            particleGun->SetOrigMomentumAmp(
                    G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
    } while ( false );
}


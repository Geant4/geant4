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
 *       Filename:  CexmcPrimaryGeneratorActionMessenger.cc
 *
 *    Description:  user assigned gun parameters
 *
 *        Version:  1.0
 *        Created:  02.11.2009 13:32:33
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include "CexmcPrimaryGeneratorAction.hh"
#include "CexmcPrimaryGeneratorActionMessenger.hh"
#include "CexmcMessenger.hh"


CexmcPrimaryGeneratorActionMessenger::CexmcPrimaryGeneratorActionMessenger(
        CexmcPrimaryGeneratorAction *  primaryGeneratorAction_ ) :
    primaryGeneratorAction( primaryGeneratorAction_ ),
    fwhmPosX( NULL ), fwhmPosY( NULL ), fwhmDirX( NULL ), fwhmDirY( NULL ),
    fwhmMomentumAmp( NULL )
{
    fwhmPosX = new G4UIcmdWithADoubleAndUnit(
                ( CexmcMessenger::gunDirName + "fwhmPosX" ).c_str(), this );
    fwhmPosX->SetGuidance( "Set positional fwhm of the beam along X axis" );
    fwhmPosX->SetParameterName( "FwhmPosX", false );
    fwhmPosX->SetDefaultValue( 0 );
    fwhmPosX->SetUnitCandidates( "mm cm m" );
    fwhmPosX->SetDefaultUnit( "cm" );
    fwhmPosX->AvailableForStates( G4State_PreInit, G4State_Idle );

    fwhmPosY = new G4UIcmdWithADoubleAndUnit(
                ( CexmcMessenger::gunDirName + "fwhmPosY" ).c_str(), this );
    fwhmPosY->SetGuidance( "Set positional fwhm of the beam along Y axis" );
    fwhmPosY->SetParameterName( "FwhmPosY", false );
    fwhmPosY->SetDefaultValue( 0 );
    fwhmPosY->SetUnitCandidates( "mm cm m" );
    fwhmPosY->SetDefaultUnit( "cm" );
    fwhmPosY->AvailableForStates( G4State_PreInit, G4State_Idle );

    fwhmDirX = new G4UIcmdWithADoubleAndUnit(
                ( CexmcMessenger::gunDirName + "fwhmDirX" ).c_str(), this );
    fwhmDirX->SetGuidance( "Set directional fwhm of the beam along X axis" );
    fwhmDirX->SetParameterName( "FwhmDirX", false );
    fwhmDirX->SetDefaultValue( 0 );
    fwhmDirX->SetUnitCandidates( "deg rad" );
    fwhmDirX->SetDefaultUnit( "deg" );
    fwhmDirX->AvailableForStates( G4State_PreInit, G4State_Idle );

    fwhmDirY = new G4UIcmdWithADoubleAndUnit(
                ( CexmcMessenger::gunDirName + "fwhmDirY" ).c_str(), this );
    fwhmDirY->SetGuidance( "Set directional fwhm of the beam along Y axis" );
    fwhmDirY->SetParameterName( "FwhmDirY", false );
    fwhmDirY->SetDefaultValue( 0 );
    fwhmDirY->SetUnitCandidates( "deg rad" );
    fwhmDirY->SetDefaultUnit( "deg" );
    fwhmDirY->AvailableForStates( G4State_PreInit, G4State_Idle );

    fwhmMomentumAmp = new G4UIcmdWithADouble(
            ( CexmcMessenger::gunDirName + "fwhmMomentumAmp" ).c_str(), this );
    fwhmMomentumAmp->SetGuidance( "Set fwhm of the beam momentum as fraction "
                                  "of its value" );
    fwhmMomentumAmp->SetParameterName( "FwhmMomentumAmp", false );
    fwhmMomentumAmp->SetDefaultValue( 0 );
    fwhmMomentumAmp->AvailableForStates( G4State_PreInit, G4State_Idle );
}


CexmcPrimaryGeneratorActionMessenger::~CexmcPrimaryGeneratorActionMessenger()
{
    delete fwhmPosX;
    delete fwhmPosY;
    delete fwhmDirX;
    delete fwhmDirY;
    delete fwhmMomentumAmp;
}


void  CexmcPrimaryGeneratorActionMessenger::SetNewValue( G4UIcommand *  cmd,
                                                         G4String  value )
{
    do
    {
        if ( cmd == fwhmPosX )
        {
            primaryGeneratorAction->SetFwhmPosX(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == fwhmPosY )
        {
            primaryGeneratorAction->SetFwhmPosY(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == fwhmDirX )
        {
            primaryGeneratorAction->SetFwhmDirX(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == fwhmDirY )
        {
            primaryGeneratorAction->SetFwhmDirY(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == fwhmMomentumAmp )
        {
            primaryGeneratorAction->SetFwhmMomentumAmp(
                        G4UIcmdWithADouble::GetNewDoubleValue( value ) );
            break;
        }
    } while ( false );
}


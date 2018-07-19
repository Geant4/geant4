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
 * =============================================================================
 *
 *       Filename:  CexmcScenePrimitivesMessenger.cc
 *
 *    Description:  draw auxiliary scene primitives
 *
 *        Version:  1.0
 *        Created:  03.01.2011 12:42:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <sstream>
#include <cctype>
#include <G4UIparameter.hh>
#include <G4UIcommand.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4Colour.hh>
#include "CexmcScenePrimitives.hh"
#include "CexmcScenePrimitivesMessenger.hh"
#include "CexmcMessenger.hh"


CexmcScenePrimitivesMessenger::CexmcScenePrimitivesMessenger(
                                    CexmcScenePrimitives *  scenePrimitives_ ) :
    scenePrimitives( scenePrimitives_ ), drawRadialLine( NULL ),
    clearRadialLines( NULL ), markTargetCenter( NULL ),
    highlightInnerCrystals( NULL ), setColour( NULL )
{
    drawRadialLine = new G4UIcmdWith3Vector(
        ( CexmcMessenger::visDirName + "drawRadialLine" ).c_str(), this );
    drawRadialLine->SetGuidance( "Draw radial line with specified theta, phi "
                                 "(both in deg!)\n    and length (in cm!)" );
    drawRadialLine->SetParameterName( "RadialLineTheta", "RadialLinePhi",
                                      "RadialLineLength", true );
    drawRadialLine->SetRange( "RadialLineLength >= 0." );
    drawRadialLine->SetDefaultValue( G4ThreeVector( 0., 0., 100. ) );
    drawRadialLine->AvailableForStates( G4State_PreInit, G4State_Idle );

    clearRadialLines = new G4UIcmdWithoutParameter(
        ( CexmcMessenger::visDirName + "clearRadialLines" ).c_str(), this );
    clearRadialLines->SetGuidance( "Clear all existing radial lines" );
    clearRadialLines->AvailableForStates( G4State_PreInit, G4State_Idle );

    markTargetCenter = new G4UIcmdWithABool(
        ( CexmcMessenger::visDirName + "markTargetCenter" ).c_str(), this );
    markTargetCenter->SetGuidance( "Mark/unmark target center" );
    markTargetCenter->SetParameterName( "MarkTargetCenter", true );
    markTargetCenter->SetDefaultValue( true );
    markTargetCenter->AvailableForStates( G4State_PreInit, G4State_Idle );

    highlightInnerCrystals = new G4UIcmdWithABool(
        ( CexmcMessenger::visDirName + "hlIC" ).c_str(), this );
    highlightInnerCrystals->SetGuidance( "Highlight inner crystals in "
                                         "calorimeters on/off" );
    highlightInnerCrystals->SetParameterName( "HighlightInnerCrystals", true );
    highlightInnerCrystals->SetDefaultValue( true );
    highlightInnerCrystals->AvailableForStates( G4State_PreInit, G4State_Idle );

    setColour = new G4UIcommand(
        ( CexmcMessenger::visDirName + "setColour" ).c_str(), this );
    setColour->SetGuidance( "Set colour of specified scene primitive" );
    G4UIparameter *  parameter( new G4UIparameter( "ScenePrimitive", 's',
                                                   false ) );
    parameter->SetGuidance( "Scene primitive, possible values:\n"
        "    tc - target center mark,\n"
        "    rl - radial lines,\n"
        "    ic - inner crystal highlights" );
    parameter->SetParameterCandidates( "tc rl ic" );
    setColour->SetParameter( parameter );
    parameter = new G4UIparameter( "red", 's', true );
    parameter->SetGuidance( "Red component or string, e.g. \"blue\", in which "
        "case succeeding colour\n    components are ignored" );
    parameter->SetDefaultValue( "1.0" );
    setColour->SetParameter( parameter );
    parameter = new G4UIparameter( "green", 'd', true );
    parameter->SetGuidance( "Green component" );
    parameter->SetDefaultValue( 1.0 );
    setColour->SetParameter( parameter );
    parameter = new G4UIparameter( "blue", 'd', true );
    parameter->SetGuidance( "Blue component" );
    parameter->SetDefaultValue( 1.0 );
    setColour->SetParameter( parameter );
    parameter = new G4UIparameter( "opacity", 'd', true );
    parameter->SetGuidance( "Opacity" );
    parameter->SetDefaultValue( 1.0 );
    setColour->SetParameter( parameter );
    setColour->AvailableForStates( G4State_PreInit, G4State_Idle );
}


CexmcScenePrimitivesMessenger::~CexmcScenePrimitivesMessenger()
{
    delete drawRadialLine;
    delete clearRadialLines;
    delete markTargetCenter;
    delete highlightInnerCrystals;
    delete setColour;
}


void  CexmcScenePrimitivesMessenger::SetNewValue( G4UIcommand *  cmd,
                                                  G4String  value )
{
    do
    {
        if ( cmd == drawRadialLine )
        {
            G4ThreeVector  line( G4UIcmdWith3Vector::GetNew3VectorValue(
                                                                    value ) );
            scenePrimitives->DrawRadialLine( line );
            break;
        }
        if ( cmd == clearRadialLines )
        {
            scenePrimitives->ClearRadialLines();
            break;
        }
        if ( cmd == markTargetCenter )
        {
            scenePrimitives->MarkTargetCenter(
                                G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
        if ( cmd == highlightInnerCrystals )
        {
            scenePrimitives->HighlightInnerCrystals(
                                G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
        if ( cmd == setColour )
        {
            G4String            name;
            G4String            redOrString;
            G4double            green( 1 );
            G4double            blue( 1 );
            G4double            opacity( 1 );
            G4Colour            colour( 1, green, blue, opacity );
            std::istringstream  iss( value );

            iss >> name >> redOrString >> green >> blue >> opacity;

            if ( std::isalpha( redOrString[ size_t( 0 ) ] ) )
            {
                G4Colour::GetColour( redOrString, colour );
            }
            else
            {
                colour = G4Colour( G4UIcommand::ConvertToDouble( redOrString ),
                                   green, blue );
            }
            colour = G4Colour( colour.GetRed(), colour.GetGreen(),
                               colour.GetBlue(), opacity );

            CexmcSPType  primitive( CexmcTargetCenterMark_SP );
            do
            {
                if ( name == "tc" )
                {
                    primitive = CexmcTargetCenterMark_SP;
                    break;
                }
                if ( name == "rl" )
                {
                    primitive = CexmcRadialLine_SP;
                    break;
                }
                if ( name == "ic" )
                {
                    primitive = CexmcInnerCrystalsHl_SP;
                    break;
                }
                return;
            } while ( false );

            scenePrimitives->SetColour( primitive, colour );
            break;
        }
    } while ( false );
}


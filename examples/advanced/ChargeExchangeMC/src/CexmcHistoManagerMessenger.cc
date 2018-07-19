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
 *       Filename:  CexmcHistoManagerMessenger.cc
 *
 *    Description:  commands to list and show histograms
 *
 *        Version:  1.0
 *        Created:  17.12.2009 21:39:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifdef CEXMC_USE_ROOT

#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithAString.hh>
#include <G4String.hh>
#include "CexmcHistoManagerMessenger.hh"
#include "CexmcHistoManager.hh"
#include "CexmcMessenger.hh"


CexmcHistoManagerMessenger::CexmcHistoManagerMessenger(
                                        CexmcHistoManager *  histoManager_ ) :
    histoManager( histoManager_ ), setVerboseLevel( NULL ), listHistos( NULL ),
    printHisto( NULL )
#ifdef CEXMC_USE_ROOTQT
    , drawHisto( NULL ), addHistoMenu( NULL ), drawHistoOptions1D( NULL ),
    drawHistoOptions2D( NULL ), drawHistoOptions3D( NULL )
#endif
{
    setVerboseLevel = new G4UIcmdWithAnInteger(
        ( CexmcMessenger::histoDirName + "verbose" ).c_str(), this );
    setVerboseLevel->SetGuidance( "0 - basic set of histograms created, "
                                  "1 - extra histograms created" );
    setVerboseLevel->SetParameterName( "Verbose", true );
    setVerboseLevel->SetDefaultValue( 0 );
    setVerboseLevel->AvailableForStates( G4State_PreInit );

    listHistos = new G4UIcmdWithoutParameter(
        ( CexmcMessenger::histoDirName + "list" ).c_str(), this );
    listHistos->SetGuidance( "List available histograms" );
    listHistos->AvailableForStates( G4State_PreInit, G4State_Idle );

    printHisto = new G4UIcmdWithAString(
        ( CexmcMessenger::histoDirName + "print" ).c_str(), this );
    printHisto->SetGuidance( "Print specified histogram" );
    printHisto->SetParameterName( "PrintHisto", false );
    printHisto->AvailableForStates( G4State_Idle );

#ifdef CEXMC_USE_ROOTQT
    drawHisto = new G4UIcmdWithAString(
        ( CexmcMessenger::histoDirName + "draw" ).c_str(), this );
    drawHisto->SetGuidance( "Draw specified histogram. The first parameter is\n"
                            "    the histogram name, the rest is draw "
                            "options.\n    Available only if the program was "
                            "launched\n    in graphical mode" );
    drawHisto->SetParameterName( "DrawHisto", false );
    drawHisto->AvailableForStates( G4State_Idle );

    addHistoMenu = new G4UIcmdWithAString(
        ( CexmcMessenger::histoDirName + "addHistoMenu" ).c_str(), this );
    addHistoMenu->SetGuidance( "Add histogram menu in GUI menu bar. The first "
                               "parameter is\n    the menu handle, the rest is "
                               "the menu label.\n    The menu cannot be added "
                               "or disabled in runtime" );
    addHistoMenu->SetParameterName( "AddHistoMenu", true );
    addHistoMenu->AvailableForStates( G4State_Idle );

    drawHistoOptions1D = new G4UIcmdWithAString(
        ( CexmcMessenger::histoDirName + "drawOptions1D" ).c_str(), this );
    drawHistoOptions1D->SetGuidance( "Set draw options for 1D histograms drawn "
                                     "from the menu.\n    The options cannot "
                                     "be changed after the menu was built");
    drawHistoOptions1D->SetParameterName( "DrawOptions1D", false );
    drawHistoOptions1D->AvailableForStates( G4State_Idle );

    drawHistoOptions2D = new G4UIcmdWithAString(
        ( CexmcMessenger::histoDirName + "drawOptions2D" ).c_str(), this );
    drawHistoOptions2D->SetGuidance( "Set draw options for 2D histograms drawn "
                                     "from the menu.\n    The options cannot "
                                     "be changed after the menu was built");
    drawHistoOptions2D->SetParameterName( "DrawOptions2D", false );
    drawHistoOptions2D->AvailableForStates( G4State_Idle );

    drawHistoOptions3D = new G4UIcmdWithAString(
        ( CexmcMessenger::histoDirName + "drawOptions3D" ).c_str(), this );
    drawHistoOptions3D->SetGuidance( "Set draw options for 3D histograms drawn "
                                     "from the menu.\n    The options cannot "
                                     "be changed after the menu was built");
    drawHistoOptions3D->SetParameterName( "DrawOptions3D", false );
    drawHistoOptions3D->AvailableForStates( G4State_Idle );
#endif
}


CexmcHistoManagerMessenger::~CexmcHistoManagerMessenger()
{
    delete setVerboseLevel;
    delete listHistos;
    delete printHisto;
#ifdef CEXMC_USE_ROOTQT
    delete drawHisto;
    delete addHistoMenu;
    delete drawHistoOptions1D;
    delete drawHistoOptions2D;
    delete drawHistoOptions3D;
#endif
}


void  CexmcHistoManagerMessenger::SetNewValue( G4UIcommand *  cmd,
                                               G4String  value )
{
    do
    {
        if ( cmd == setVerboseLevel )
        {
            histoManager->SetVerboseLevel(
                                G4UIcmdWithAnInteger::GetNewIntValue( value ) );
            break;
        }
        if ( cmd == listHistos )
        {
            histoManager->List();
            break;
        }
        if ( cmd == printHisto )
        {
            histoManager->Print( value );
            break;
        }
#ifdef CEXMC_USE_ROOTQT
        if ( cmd == drawHisto )
        {
            size_t  delimPos( value.find_first_of( " \t" ) );
            size_t  delimPosEnd( G4String::npos );
            if ( delimPos != G4String::npos )
                delimPosEnd = value.find_first_not_of( " \t", delimPos );
            histoManager->Draw( std::string( value, 0, delimPos ),
                                delimPosEnd == G4String::npos ? "" :
                                                value.c_str() + delimPosEnd );
            break;
        }
        if ( cmd == addHistoMenu )
        {
            size_t  delimPos( value.find_first_of( " \t" ) );
            size_t  delimPosEnd( G4String::npos );
            if ( delimPos != G4String::npos )
                delimPosEnd = value.find_first_not_of( " \t", delimPos );
            histoManager->AddHistoMenu( std::string( value, 0, delimPos ),
                                        delimPosEnd == G4String::npos ? "" :
                                                value.c_str() + delimPosEnd );
            break;
        }
        if ( cmd == drawHistoOptions1D )
        {
            histoManager->SetDrawOptions1D( value );
            break;
        }
        if ( cmd == drawHistoOptions2D )
        {
            histoManager->SetDrawOptions2D( value );
            break;
        }
        if ( cmd == drawHistoOptions3D )
        {
            histoManager->SetDrawOptions3D( value );
            break;
        }
#endif
    } while ( false );
}

#endif


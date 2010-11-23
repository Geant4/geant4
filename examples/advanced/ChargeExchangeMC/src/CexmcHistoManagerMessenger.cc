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
#include <G4UIcmdWithAString.hh>
#include <G4String.hh>
#include "CexmcHistoManagerMessenger.hh"
#include "CexmcHistoManager.hh"
#include "CexmcMessenger.hh"


CexmcHistoManagerMessenger::CexmcHistoManagerMessenger() :
    listHistos( NULL ), printHisto( NULL )
#ifdef CEXMC_USE_ROOTQT
    , drawHisto( NULL )
#endif
{
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
                            "    the histogram name, the second is draw "
                            "options.\n    Available only if the program was "
                            "launched\n    in graphical mode" );
    drawHisto->SetParameterName( "DrawHisto", false );
    drawHisto->AvailableForStates( G4State_Idle );
#endif
}


CexmcHistoManagerMessenger::~CexmcHistoManagerMessenger()
{
    delete listHistos;
    delete printHisto;
#ifdef CEXMC_USE_ROOTQT
    delete drawHisto;
#endif
}


void CexmcHistoManagerMessenger::SetNewValue(
                                        G4UIcommand *  cmd, G4String  value )
{
    CexmcHistoManager *  histoManager( CexmcHistoManager::Instance() );

    do
    {
        if ( cmd == listHistos )
        {
            histoManager->List();
        }
        if ( cmd == printHisto )
        {
            histoManager->Print( value );
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
        }
#endif
    } while ( false );
}

#endif


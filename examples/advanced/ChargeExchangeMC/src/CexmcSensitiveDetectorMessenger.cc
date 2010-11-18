/*
 * ============================================================================
 *
 *       Filename:  CexmcSensitiveDetectorMessenger.cc
 *
 *    Description:  sensitive detector messenger (verbose level etc.)
 *
 *        Version:  1.0
 *        Created:  15.11.2009 14:10:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4UIdirectory.hh>
#include <G4UImanager.hh>
#include <G4UIcommandTree.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4VPrimitiveScorer.hh>
#include "CexmcSensitiveDetectorMessenger.hh"
#include "CexmcMessenger.hh"


CexmcSensitiveDetectorMessenger::CexmcSensitiveDetectorMessenger(
            G4VPrimitiveScorer *  scorer, const G4String &  detectorName ) :
    scorer( scorer ), detectorPath( NULL ), setVerboseLevel( NULL )
{
    G4String       detectorFullPath(
            ( CexmcMessenger::detectorDirName + detectorName + "/" ).c_str() );
    G4UImanager *  uiManager( G4UImanager::GetUIpointer() );
    if ( ! uiManager->GetTree()->FindCommandTree( detectorFullPath.c_str() ) )
    {
        detectorPath = new G4UIdirectory( detectorFullPath.c_str() );
        detectorPath->SetGuidance( "Settings for given sensitive detector" );
    }
    setVerboseLevel = new G4UIcmdWithAnInteger(
                            ( detectorFullPath + "verbose" ).c_str(), this );
    setVerboseLevel->SetGuidance( "0 - do not print messages, "
                                  "1 - print messages after an event" );
    setVerboseLevel->SetParameterName( "Verbose", true );
    setVerboseLevel->SetDefaultValue( 0 );
    setVerboseLevel->AvailableForStates( G4State_PreInit, G4State_Idle );

}


CexmcSensitiveDetectorMessenger::~CexmcSensitiveDetectorMessenger()
{
    delete detectorPath;
    delete setVerboseLevel;
}


void  CexmcSensitiveDetectorMessenger::SetNewValue( G4UIcommand *  cmd,
                                                    G4String  value )
{
    do
    {
        if ( cmd == setVerboseLevel )
        {
            scorer->SetVerboseLevel(
                                G4UIcmdWithAnInteger::GetNewIntValue( value ) );
            break;
        }
    } while ( false );
}


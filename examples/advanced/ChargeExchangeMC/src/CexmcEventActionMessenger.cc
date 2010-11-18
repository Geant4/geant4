/*
 * ============================================================================
 *
 *       Filename:  CexmcEventActionMessenger.cc
 *
 *    Description:  event action messenger (verbose level etc.)
 *
 *        Version:  1.0
 *        Created:  25.11.2009 14:45:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithABool.hh>
#include "CexmcEventActionMessenger.hh"
#include "CexmcEventAction.hh"
#include "CexmcMessenger.hh"


CexmcEventActionMessenger::CexmcEventActionMessenger(
                                            CexmcEventAction *  eventAction ) :
    eventAction( eventAction ), setVerboseLevel( NULL ),
    setVerboseDrawLevel( NULL ), drawTrajectoryMarkers( NULL )
{
    setVerboseLevel = new G4UIcmdWithAnInteger(
           ( CexmcMessenger::eventDirName + "verbose" ).c_str() , this );
    setVerboseLevel->SetGuidance( "\n    0 - do not print messages,\n"
                        "    1 - print messages if studied interaction "
                                "triggered,\n"
                        "    2 - print messages on event trigger\n"
                        "        (reconstructed data will also be printed in "
                                "this case),\n"
                        "    3 - print messages if studied interaction or "
                                "event\n        triggered,\n"
                        "    4 - print messages on every event" );
    setVerboseLevel->SetParameterName( "Verbose", false );
    setVerboseLevel->SetDefaultValue( 0 );
    setVerboseLevel->SetRange( "Verbose >= 0 && Verbose <= 4" );
    setVerboseLevel->AvailableForStates( G4State_PreInit, G4State_Idle );

    setVerboseDrawLevel = new G4UIcmdWithAnInteger(
           ( CexmcMessenger::visDirName + "verbose" ).c_str() , this );
    setVerboseDrawLevel->SetGuidance( "\n    0 - draw nothing,\n"
                        "    1 - draw trajectories and track points if studied "
                                "interaction\n        triggered,\n"
                        "    2 - draw trajectories and track points on event "
                                "trigger\n"
                        "        (reconstructed data will also be drawn in "
                                "this case),\n"
                        "    3 - draw trajectories and track points if "
                                "studied interaction\n        or event "
                                "triggered,\n"
                        "    4 - draw trajectories and/or track points on "
                                "every event" );
    setVerboseDrawLevel->SetParameterName( "VerboseDraw", false );
    setVerboseDrawLevel->SetDefaultValue( 4 );
    setVerboseDrawLevel->SetRange( "VerboseDraw >= 0 && VerboseDraw <= 4" );
    setVerboseDrawLevel->AvailableForStates( G4State_PreInit, G4State_Idle );

    drawTrajectoryMarkers = new G4UIcmdWithABool(
           ( CexmcMessenger::visDirName + "drawTrajMarkers" ).c_str() , this );
    drawTrajectoryMarkers->SetGuidance( "draw markers in places where "
                                        "trajectories change" );
    drawTrajectoryMarkers->SetParameterName( "DrawTrajMarkers", false );
    drawTrajectoryMarkers->SetDefaultValue( false );
    drawTrajectoryMarkers->AvailableForStates( G4State_PreInit, G4State_Idle );
}


CexmcEventActionMessenger::~CexmcEventActionMessenger()
{
    delete setVerboseLevel;
    delete setVerboseDrawLevel;
    delete drawTrajectoryMarkers;
}


void  CexmcEventActionMessenger::SetNewValue( G4UIcommand *  cmd,
                                              G4String  value )
{
    do
    {
        if ( cmd == setVerboseLevel )
        {
            eventAction->SetVerboseOnCexmcLevel(
                                G4UIcmdWithAnInteger::GetNewIntValue( value ) );
            break;
        }
        if ( cmd == setVerboseDrawLevel )
        {
            eventAction->SetVerboseDrawLevel(
                                G4UIcmdWithAnInteger::GetNewIntValue( value ) );
            break;
        }
        if ( cmd == drawTrajectoryMarkers )
        {
            eventAction->DrawTrajectoryMarkers(
                                G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
    } while ( false );
}


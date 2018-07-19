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
#include "CexmcEventActionMessenger.hh"
#include "CexmcEventAction.hh"
#include "CexmcMessenger.hh"


CexmcEventActionMessenger::CexmcEventActionMessenger(
                                            CexmcEventAction *  eventAction_ ) :
    eventAction( eventAction_ ), setVerboseLevel( NULL ),
    setVerboseDrawLevel( NULL )
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
    setVerboseLevel->SetRange( "Verbose >= 0 && Verbose <= 4" );
    setVerboseLevel->SetDefaultValue( 0 );
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
    setVerboseDrawLevel->SetRange( "VerboseDraw >= 0 && VerboseDraw <= 4" );
    setVerboseDrawLevel->SetDefaultValue( 4 );
    setVerboseDrawLevel->AvailableForStates( G4State_PreInit, G4State_Idle );
}


CexmcEventActionMessenger::~CexmcEventActionMessenger()
{
    delete setVerboseLevel;
    delete setVerboseDrawLevel;
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
    } while ( false );
}


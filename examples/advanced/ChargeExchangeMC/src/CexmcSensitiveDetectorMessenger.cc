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
#include <G4MultiFunctionalDetector.hh>
#include "CexmcSensitiveDetectorMessenger.hh"
#include "CexmcMessenger.hh"
#include "CexmcException.hh"


CexmcSensitiveDetectorMessenger::CexmcSensitiveDetectorMessenger(
                                            G4VPrimitiveScorer *  scorer_ ) :
    scorer( scorer_ ), detectorPath( NULL ), setVerboseLevel( NULL )
{
    G4MultiFunctionalDetector *  detector(
                                        scorer->GetMultiFunctionalDetector() );
    /* detector of the scorer must have been already initialized provided
     * CexmcPrimitiveScorer::InitializeMessenger() was properly called upon
     * creation of the setup */
    if ( ! detector )
        throw CexmcException( CexmcWeirdException );

    G4String       detectorFullPath( ( CexmcMessenger::detectorDirName +
                       detector->GetName() + "/" + scorer->GetName() + "/" ).
                                                                    c_str() );
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


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
 *       Filename:  CexmcPhysicsManagerMessenger.cc
 *
 *    Description:  physics manager messenger (max IL correction etc.)
 *
 *        Version:  1.0
 *        Created:  16.10.2010 14:15:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <G4UIcmdWithADoubleAndUnit.hh>
#include "CexmcPhysicsManager.hh"
#include "CexmcPhysicsManagerMessenger.hh"
#include "CexmcMessenger.hh"


CexmcPhysicsManagerMessenger::CexmcPhysicsManagerMessenger(
                                    CexmcPhysicsManager *  physicsManager_ ) :
    physicsManager( physicsManager_ ), setMaxILCorrection( NULL )
{
    setMaxILCorrection = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::physicsDirName + "setMaxILCorrection" ).c_str(),
        this );
    setMaxILCorrection->SetGuidance( "Correction of maximum interaction "
        "length in the\n    target; can be used for precise tuning of the "
        "center of\n    distribution of interaction vertices in the target" );
    setMaxILCorrection->SetParameterName( "MaxILCorrection", false );
    setMaxILCorrection->SetUnitCandidates( "mm cm m" );
    setMaxILCorrection->SetDefaultUnit( "mm" );
    setMaxILCorrection->AvailableForStates( G4State_PreInit, G4State_Idle );
}


CexmcPhysicsManagerMessenger::~CexmcPhysicsManagerMessenger()
{
    delete setMaxILCorrection;
}


void  CexmcPhysicsManagerMessenger::SetNewValue( G4UIcommand *  cmd,
                                                 G4String  value )
{
    do
    {
        if ( cmd == setMaxILCorrection )
        {
            physicsManager->SetMaxILCorrection(
                    G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
    } while ( false );
}


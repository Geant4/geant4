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
                                    CexmcPhysicsManager *  physicsManager ) :
    physicsManager( physicsManager ), setMaxILCorrection( NULL )
{
    setMaxILCorrection = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::physicsDirName + "setMaxILCorrection" ).c_str(),
        this );
    setMaxILCorrection->SetGuidance( "Correction of maximum interaction "
        "length in the\n    target; can be used for precise tuning of the "
        "center of\n    distribution of interaction vertices in the target" );
    setMaxILCorrection->SetParameterName( "MaxILCorrection", false );
    setMaxILCorrection->SetDefaultUnit( "mm" );
    setMaxILCorrection->SetUnitCandidates( "mm cm m" );
    setMaxILCorrection->AvailableForStates( G4State_PreInit, G4State_Idle );
}


CexmcPhysicsManagerMessenger::~CexmcPhysicsManagerMessenger()
{
    delete setMaxILCorrection;
}


void CexmcPhysicsManagerMessenger::SetNewValue( G4UIcommand *  cmd,
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


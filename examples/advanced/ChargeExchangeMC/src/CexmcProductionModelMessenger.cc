/*
 * ============================================================================
 *
 *       Filename:  CexmcProductionModelMessenger.cc
 *
 *    Description:  set various production model aspects
 *
 *        Version:  1.0
 *        Created:  03.11.2009 16:01:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWith3Vector.hh>
#include "CexmcProductionModel.hh"
#include "CexmcProductionModelMessenger.hh"
#include "CexmcMessenger.hh"


CexmcProductionModelMessenger::CexmcProductionModelMessenger(
                                    CexmcProductionModel *  productionModel ) :
    productionModel( productionModel ), applyFermiMotion( NULL ),
    setAngularRange( NULL ), addAngularRange( NULL )
{
    applyFermiMotion = new G4UIcmdWithABool(
        ( CexmcMessenger::physicsDirName + "applyFermiMotionInTarget" ).c_str(),
        this );
    applyFermiMotion->SetGuidance( "Switch on/off fermi motion in target "
                                   "nuclei" );
    applyFermiMotion->SetParameterName( "ApplyFermiMotionInTarget", false );
    applyFermiMotion->SetDefaultValue( false );
    applyFermiMotion->AvailableForStates( G4State_PreInit, G4State_Idle );

    setAngularRange = new G4UIcmdWith3Vector(
        ( CexmcMessenger::physicsDirName + "setAngularRange" ).c_str(), this );
    setAngularRange->SetGuidance(
        "\n    Set angular range of interest given in values of cosinus;\n"
        "    first two values give the range (descending or ascending),\n"
        "    third value gives number of equal divisions within the\n"
        "    range." );
    setAngularRange->SetParameterName( "ARangeTop", "ARangeBottom",
                                       "ARangeNmbOfDivisions", false );
    setAngularRange->SetRange(
        "ARangeTop >= -1.0 && ARangeTop <= 1.0 && "
        "ARangeBottom >= -1.0 && ARangeBottom <= 1.0 && "
        "ARangeNmbOfDivisions >= 1" );
    setAngularRange->AvailableForStates( G4State_PreInit, G4State_Idle );

    addAngularRange = new G4UIcmdWith3Vector(
        ( CexmcMessenger::physicsDirName + "addAngularRange" ).c_str(), this );
    addAngularRange->SetGuidance(
        "\n    Add angular range of interest given in values of cosinus;\n"
        "    first two values give the range (descending or ascending),\n"
        "    third value gives number of equal divisions within the\n"
        "    range." );
    addAngularRange->SetParameterName( "ARangeTop", "ARangeBottom",
                                       "ARangeNmbOfDivisions", false );
    addAngularRange->SetRange(
        "ARangeTop >= -1.0 && ARangeTop <= 1.0 && "
        "ARangeBottom >= -1.0 && ARangeBottom <= 1.0 && "
        "ARangeNmbOfDivisions >= 1" );
    addAngularRange->AvailableForStates( G4State_PreInit, G4State_Idle );
}


CexmcProductionModelMessenger::~CexmcProductionModelMessenger()
{
    delete applyFermiMotion;
    delete setAngularRange;
    delete addAngularRange;
}


void CexmcProductionModelMessenger::SetNewValue( G4UIcommand *  cmd,
                                                 G4String  value )
{
    do
    {
        if ( cmd == applyFermiMotion )
        {
            productionModel->ApplyFermiMotion(
                                G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
        if ( cmd == setAngularRange )
        {
            G4ThreeVector  vec( G4UIcmdWith3Vector::GetNew3VectorValue(
                                                                    value ) );
            G4double       top( std::max( vec.x(), vec.y() ) );
            G4double       bottom( std::min( vec.x(), vec.y() ) );
            productionModel->SetAngularRange( top, bottom,
                                              static_cast< int >( vec.z() ) );
            break;
        }
        if ( cmd == addAngularRange )
        {
            G4ThreeVector  vec( G4UIcmdWith3Vector::GetNew3VectorValue(
                                                                    value ) );
            G4double       top( std::max( vec.x(), vec.y() ) );
            G4double       bottom( std::min( vec.x(), vec.y() ) );
            productionModel->AddAngularRange( top, bottom,
                                              static_cast< int >( vec.z() ) );
            break;
        }
    } while ( false );
}


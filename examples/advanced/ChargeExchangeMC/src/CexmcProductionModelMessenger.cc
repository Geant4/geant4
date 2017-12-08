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
                                    CexmcProductionModel *  productionModel_ ) :
    productionModel( productionModel_ ), applyFermiMotion( NULL ),
    setAngularRange( NULL ), addAngularRange( NULL )
{
    applyFermiMotion = new G4UIcmdWithABool(
        ( CexmcMessenger::physicsDirName + "applyFermiMotionInTarget" ).c_str(),
        this );
    applyFermiMotion->SetGuidance( "Switch on/off fermi motion in target "
                                   "nuclei" );
    applyFermiMotion->SetParameterName( "ApplyFermiMotionInTarget", true );
    applyFermiMotion->SetDefaultValue( true );
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


void  CexmcProductionModelMessenger::SetNewValue( G4UIcommand *  cmd,
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


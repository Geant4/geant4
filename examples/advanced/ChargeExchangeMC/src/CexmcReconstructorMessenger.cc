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
 *       Filename:  CexmcReconstructorMessenger.cc
 *
 *    Description:  reconstructor messenger
 *
 *        Version:  1.0
 *        Created:  02.12.2009 15:38:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include "CexmcReconstructorMessenger.hh"
#include "CexmcReconstructor.hh"
#include "CexmcMessenger.hh"
#include "CexmcCommon.hh"


CexmcReconstructorMessenger::CexmcReconstructorMessenger(
                                        CexmcReconstructor *  reconstructor_ ) :
    reconstructor( reconstructor_ ),
    setCalorimeterEntryPointDefinitionAlgorithm( NULL ),
    setCrystalSelectionAlgorithm( NULL ), useInnerRefCrystal( NULL ),
    setCalorimeterEntryPointDepth( NULL )
{
    setCalorimeterEntryPointDefinitionAlgorithm = new G4UIcmdWithAString(
        ( CexmcMessenger::reconstructorDirName + "entryPointDefinitionAlgo" ).
                c_str(), this );
    setCalorimeterEntryPointDefinitionAlgorithm->SetGuidance(
        "\n    Algorithm to reconstruct entry point of output particle"
        "\n    decay products in calorimeter"
        "\n    (none of the following algorithms reconstruct directions)\n"
        "    center - entry points defined in the center of the\n"
        "             calorimeters,\n"
        "    simple - entry points defined in the center of the crystal\n"
        "             that has maximum energy deposit value,\n"
        "    linear - entry points defined by linear weights of energy\n"
        "             deposit in crystals,\n"
        "    sqrt - entry points defined by square root weights of\n"
        "           energy deposit in crystals" );
    setCalorimeterEntryPointDefinitionAlgorithm->SetParameterName(
                                            "EntryPointDefinitionAlgo", false );
    setCalorimeterEntryPointDefinitionAlgorithm->SetCandidates(
                                            "center simple linear sqrt" );
    setCalorimeterEntryPointDefinitionAlgorithm->SetDefaultValue( "sqrt" );
    setCalorimeterEntryPointDefinitionAlgorithm->AvailableForStates(
                                            G4State_PreInit, G4State_Idle );

    setCalorimeterEntryPointDepthDefinitionAlgorithm = new G4UIcmdWithAString(
        ( CexmcMessenger::reconstructorDirName +
                              "entryPointDepthDefinitionAlgo" ).c_str(), this );
    setCalorimeterEntryPointDepthDefinitionAlgorithm->SetGuidance(
        "\n    Algorithm to reconstruct entry point depth of output\n"
        "    particle decay products in calorimeter\n"
        "    (value is defined by 'entryPointDepth' parameter)\n"
        "    plain - depth is a constant\n"
        "    sphere - depth depends on X and Y of calorimeter entry\n"
        "             points and locates on surface of a sphere\n"
        "             with origin in the center of the target;\n"
        "             radius of the sphere is sum of distance to\n"
        "             the calorimeter and 'entryPointDepth' value" );
    setCalorimeterEntryPointDepthDefinitionAlgorithm->SetParameterName(
                                    "EntryPointDepthDefinitionAlgo", false );
    setCalorimeterEntryPointDepthDefinitionAlgorithm->SetCandidates(
                                    "plain sphere" );
    setCalorimeterEntryPointDepthDefinitionAlgorithm->SetDefaultValue(
                                    "plain" );
    setCalorimeterEntryPointDepthDefinitionAlgorithm->AvailableForStates(
                                                G4State_PreInit, G4State_Idle );

    setCrystalSelectionAlgorithm = new G4UIcmdWithAString(
        ( CexmcMessenger::reconstructorDirName + "crystalSelectionAlgo" ).
                c_str(), this );
    setCrystalSelectionAlgorithm->SetGuidance(
        "\n    Choose crystals to be selected in weighted entry point\n"
        "    reconstruction algorithms\n"
        "    all - all,\n"
        "    adjacent - crystal with maximum energy deposit and\n"
        "               adjacent crystals" );
    setCrystalSelectionAlgorithm->SetParameterName( "CrystalSelAlgo", false );
    setCrystalSelectionAlgorithm->SetCandidates( "all adjacent" );
    setCrystalSelectionAlgorithm->SetDefaultValue( "all" );
    setCrystalSelectionAlgorithm->AvailableForStates( G4State_PreInit,
                                                      G4State_Idle );

    useInnerRefCrystal = new G4UIcmdWithABool(
        ( CexmcMessenger::reconstructorDirName + "useInnerRefCrystal" ).
                c_str(), this );
    useInnerRefCrystal->SetGuidance(
        "\n    Defines that if the crystal with maximum energy deposit in\n"
        "    calorimeter is an outer crystal then the closest inner crystal\n"
        "    will be chosen as the reference for adjacent crystal selection\n"
        "    algorithm and simple entry point definition algorithm. It also\n"
        "    affects energy deposit collection if adjacent crystals\n"
        "    algorithm was chosen for that. If not set then the reference\n"
        "    crystal will be found from all crystals in calorimeter" );
    useInnerRefCrystal->SetParameterName( "UseInnerRefCrystal", true );
    useInnerRefCrystal->SetDefaultValue( true );
    useInnerRefCrystal->AvailableForStates( G4State_PreInit, G4State_Idle );

    setCalorimeterEntryPointDepth = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "entryPointDepth" ).c_str(),
        this );
    setCalorimeterEntryPointDepth->SetGuidance(
        "\n    Depth of entry point used in reconstruction of angle\n"
        "    between output particle decay products" );
    setCalorimeterEntryPointDepth->SetParameterName( "EntryPointDepth", false );
    setCalorimeterEntryPointDepth->SetDefaultValue( 0 );
    setCalorimeterEntryPointDepth->SetUnitCandidates( "mm cm m" );
    setCalorimeterEntryPointDepth->SetDefaultUnit( "cm" );
    setCalorimeterEntryPointDepth->AvailableForStates( G4State_PreInit,
                                                       G4State_Idle );
}


CexmcReconstructorMessenger::~CexmcReconstructorMessenger()
{
    delete setCalorimeterEntryPointDefinitionAlgorithm;
    delete setCalorimeterEntryPointDepthDefinitionAlgorithm;
    delete setCrystalSelectionAlgorithm;
    delete useInnerRefCrystal;
    delete setCalorimeterEntryPointDepth;
}


void  CexmcReconstructorMessenger::SetNewValue( G4UIcommand *  cmd,
                                                G4String  value )
{
    do
    {
        if ( cmd == setCalorimeterEntryPointDefinitionAlgorithm )
        {
            CexmcCalorimeterEntryPointDefinitionAlgorithm
                            epDefinitionAlgorithm( CexmcEntryPointInTheCenter );
            do
            {
                if ( value == "simple" )
                {
                    epDefinitionAlgorithm =
                            CexmcEntryPointInTheCenterOfCrystalWithMaxED;
                    break;
                }
                if ( value == "linear" )
                {
                    epDefinitionAlgorithm = CexmcEntryPointByLinearEDWeights;
                    break;
                }
                if ( value == "sqrt" )
                {
                    epDefinitionAlgorithm = CexmcEntryPointBySqrtEDWeights;
                    break;
                }
            } while ( false );
            reconstructor->SetCalorimeterEntryPointDefinitionAlgorithm(
                                                        epDefinitionAlgorithm );
            break;
        }
        if ( cmd == setCalorimeterEntryPointDepthDefinitionAlgorithm )
        {
            CexmcCalorimeterEntryPointDepthDefinitionAlgorithm
                        epDepthDefinitionAlgorithm( CexmcEntryPointDepthPlain );
            do
            {
                if ( value == "sphere" )
                {
                    epDepthDefinitionAlgorithm = CexmcEntryPointDepthSphere;
                    break;
                }
            } while ( false );
            reconstructor->SetCalorimeterEntryPointDepthDefinitionAlgorithm(
                                                epDepthDefinitionAlgorithm );
            break;
        }
        if ( cmd == setCrystalSelectionAlgorithm )
        {
            CexmcCrystalSelectionAlgorithm
                                        csAlgorithm( CexmcSelectAllCrystals );
            do
            {
                if ( value == "adjacent" )
                {
                    csAlgorithm = CexmcSelectAdjacentCrystals;
                    break;
                }
            } while ( false );
            reconstructor->SetCrystalSelectionAlgorithm( csAlgorithm );
            break;
        }
        if ( cmd == useInnerRefCrystal )
        {
            reconstructor->UseInnerRefCrystal(
                        G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
        if ( cmd == setCalorimeterEntryPointDepth )
        {
            reconstructor->SetCalorimeterEntryPointDepth(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
    } while ( false );
}


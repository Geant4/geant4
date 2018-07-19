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
 *       Filename:  CexmcChargeExchangeReconstructorMessenger.cc
 *
 *    Description:  charge exchange reconstructor messenger
 *
 *        Version:  1.0
 *        Created:  14.12.2009 17:53:33
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
#include "CexmcChargeExchangeReconstructorMessenger.hh"
#include "CexmcChargeExchangeReconstructor.hh"
#include "CexmcMessenger.hh"


CexmcChargeExchangeReconstructorMessenger::
        CexmcChargeExchangeReconstructorMessenger(
                        CexmcChargeExchangeReconstructor *  reconstructor_ ) :
            reconstructor( reconstructor_ ), useTableMass( NULL ),
            useMassCut( NULL ), mCutOPCenter( NULL ), mCutNOPCenter( NULL ),
            mCutOPWidth( NULL ), mCutNOPWidth( NULL ), mCutAngle( NULL ),
            useAbsorbedEnergyCut( NULL ), aeCutCLCenter( NULL ),
            aeCutCRCenter( NULL ), aeCutCLWidth( NULL ), aeCutCRWidth( NULL ),
            aeCutAngle( NULL ), setExpectedMomentumAmp( NULL ),
            setExpectedMomentumAmpDiff( NULL ), setEDCollectionAlgorithm( NULL )
{
    useTableMass = new G4UIcmdWithABool(
        ( CexmcMessenger::reconstructorDirName + "useTableMass" ).c_str(),
        this );
    useTableMass->SetGuidance( "\n    If true then reconstructor will use "
        "table mass of output\n    particle when building output particle "
        "energy,\n    otherwise reconstructed mass will be used" );
    useTableMass->SetParameterName( "UseTableMass", true );
    useTableMass->SetDefaultValue( true );
    useTableMass->AvailableForStates( G4State_PreInit, G4State_Idle );

    useMassCut = new G4UIcmdWithABool(
        ( CexmcMessenger::reconstructorDirName + "useMassCut" ).c_str(), this );
    useMassCut->SetGuidance( "\n    Use elliptical cut for masses of output "
                             "particle\n    and nucleus output particle" );
    useMassCut->SetParameterName( "UseMassCut", true );
    useMassCut->SetDefaultValue( true );
    useMassCut->AvailableForStates( G4State_PreInit, G4State_Idle );

    mCutOPCenter = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "mCutOPCenter" ).c_str(),
        this );
    mCutOPCenter->SetGuidance( "\n    Center of the ellipse in output particle "
                               "mass coordinate" );
    mCutOPCenter->SetParameterName( "MCutOPCenter", false );
    mCutOPCenter->SetDefaultValue( reconstructor->GetProductionModelData().
                                   outputParticle->GetPDGMass() );
    mCutOPCenter->SetUnitCandidates( "eV keV MeV GeV" );
    mCutOPCenter->SetDefaultUnit( "MeV" );
    mCutOPCenter->AvailableForStates( G4State_PreInit, G4State_Idle );

    mCutNOPCenter = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "mCutNOPCenter" ).c_str(),
        this );
    mCutNOPCenter->SetGuidance( "\n    Center of the ellipse in nucleus output "
                                "particle mass\n    coordinate" );
    mCutNOPCenter->SetParameterName( "MCutNOPCenter", false );
    mCutNOPCenter->SetDefaultValue( reconstructor->GetProductionModelData().
                                    nucleusOutputParticle->GetPDGMass() );
    mCutNOPCenter->SetUnitCandidates( "eV keV MeV GeV" );
    mCutNOPCenter->SetDefaultUnit( "MeV" );
    mCutNOPCenter->AvailableForStates( G4State_PreInit, G4State_Idle );

    mCutOPWidth = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "mCutOPWidth" ).c_str(),
        this );
    mCutOPWidth->SetGuidance( "\n    Width of the ellipse in output particle "
                               "mass coordinate" );
    mCutOPWidth->SetParameterName( "MCutOPWidth", false );
    mCutOPWidth->SetDefaultValue( reconstructor->GetProductionModelData().
                                  outputParticle->GetPDGMass() * 0.1 );
    mCutOPWidth->SetUnitCandidates( "eV keV MeV GeV" );
    mCutOPWidth->SetDefaultUnit( "MeV" );
    mCutOPWidth->AvailableForStates( G4State_PreInit, G4State_Idle );

    mCutNOPWidth = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "mCutNOPWidth" ).c_str(),
        this );
    mCutNOPWidth->SetGuidance( "\n    Width of the ellipse in nucleus output "
                               "particle mass\n    coordinate" );
    mCutNOPWidth->SetParameterName( "MCutNOPWidth", false );
    mCutNOPWidth->SetDefaultValue( reconstructor->GetProductionModelData().
                                   nucleusOutputParticle->GetPDGMass() * 0.1 );
    mCutNOPWidth->SetUnitCandidates( "eV keV MeV GeV" );
    mCutNOPWidth->SetDefaultUnit( "MeV" );
    mCutNOPWidth->AvailableForStates( G4State_PreInit, G4State_Idle );

    mCutAngle = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "mCutAngle" ).c_str(),
        this );
    mCutAngle->SetGuidance( "\n    Angle of the ellipse" );
    mCutAngle->SetParameterName( "MCutAngle", false );
    mCutAngle->SetDefaultValue( 0 );
    mCutAngle->SetUnitCandidates( "deg rad" );
    mCutAngle->SetDefaultUnit( "deg" );
    mCutAngle->AvailableForStates( G4State_PreInit, G4State_Idle );

    useAbsorbedEnergyCut = new G4UIcmdWithABool(
        ( CexmcMessenger::reconstructorDirName + "useAbsorbedEnergyCut" ).
            c_str(), this );
    useAbsorbedEnergyCut->SetGuidance( "\n    Use elliptical cut for absorbed "
                                       "energies in\n    calorimeters" );
    useAbsorbedEnergyCut->SetParameterName( "UseAbsorbedEnergyCut", true );
    useAbsorbedEnergyCut->SetDefaultValue( true );
    useAbsorbedEnergyCut->AvailableForStates( G4State_PreInit, G4State_Idle );

    aeCutCLCenter = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "aeCutCLCenter" ).c_str(),
        this );
    aeCutCLCenter->SetGuidance( "\n    Center of the ellipse in left "
                                "calorimeter\n    absorbed energy coordinate" );
    aeCutCLCenter->SetParameterName( "AECutCLCenter", false );
    aeCutCLCenter->SetDefaultValue( 0 );
    aeCutCLCenter->SetUnitCandidates( "eV keV MeV GeV" );
    aeCutCLCenter->SetDefaultUnit( "MeV" );
    aeCutCLCenter->AvailableForStates( G4State_PreInit, G4State_Idle );

    aeCutCRCenter = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "aeCutCRCenter" ).c_str(),
        this );
    aeCutCRCenter->SetGuidance( "\n    Center of the ellipse in right "
                                "calorimeter\n    absorbed energy coordinate" );
    aeCutCRCenter->SetParameterName( "AECutCRCenter", false );
    aeCutCRCenter->SetDefaultValue( 0 );
    aeCutCRCenter->SetUnitCandidates( "eV keV MeV GeV" );
    aeCutCRCenter->SetDefaultUnit( "MeV" );
    aeCutCRCenter->AvailableForStates( G4State_PreInit, G4State_Idle );

    aeCutCLWidth = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "aeCutCLWidth" ).c_str(),
        this );
    aeCutCLWidth->SetGuidance( "\n    Width of the ellipse in left calorimeter"
                               "\n    absorbed energy coordinate" );
    aeCutCLWidth->SetParameterName( "AECutCLWidth", false );
    aeCutCLWidth->SetDefaultValue( 0 );
    aeCutCLWidth->SetUnitCandidates( "eV keV MeV GeV" );
    aeCutCLWidth->SetDefaultUnit( "MeV" );
    aeCutCLWidth->AvailableForStates( G4State_PreInit, G4State_Idle );

    aeCutCRWidth = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "aeCutCRWidth" ).c_str(),
        this );
    aeCutCRWidth->SetGuidance( "\n    Width of the ellipse in right calorimeter"
                               "\n    absorbed energy coordinate" );
    aeCutCRWidth->SetParameterName( "AECutCRWidth", false );
    aeCutCRWidth->SetDefaultValue( 0 );
    aeCutCRWidth->SetUnitCandidates( "eV keV MeV GeV" );
    aeCutCRWidth->SetDefaultUnit( "MeV" );
    aeCutCRWidth->AvailableForStates( G4State_PreInit, G4State_Idle );

    aeCutAngle = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "aeCutAngle" ).c_str(),
        this );
    aeCutAngle->SetGuidance( "\n    Angle of the ellipse" );
    aeCutAngle->SetParameterName( "AECutAngle", false );
    aeCutAngle->SetDefaultValue( 0 );
    aeCutAngle->SetUnitCandidates( "deg rad" );
    aeCutAngle->SetDefaultUnit( "deg" );
    aeCutAngle->AvailableForStates( G4State_PreInit, G4State_Idle );

    setExpectedMomentumAmp = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "momentumAmp" ).c_str(),
        this );
    setExpectedMomentumAmp->SetGuidance( "\n    Momentum of a beam particle "
        "expected in the target;\n    this value may differ from original "
        "momentum\n    of the beam as far as profile data of the beam refer\n"
        "    to the place where it starts. This parameter is used only\n"
        "    in reconstruction procedure");
    setExpectedMomentumAmp->SetParameterName( "RecMomentumAmp", false );
    setExpectedMomentumAmp->SetRange( "RecMomentumAmp > 0" );
    setExpectedMomentumAmp->SetUnitCandidates( "eV keV MeV GeV" );
    setExpectedMomentumAmp->SetDefaultUnit( "MeV" );
    setExpectedMomentumAmp->AvailableForStates( G4State_PreInit, G4State_Idle );

    setExpectedMomentumAmpDiff = new G4UIcmdWithADoubleAndUnit(
        ( CexmcMessenger::reconstructorDirName + "momentumAmpDiff" ).c_str(),
        this );
    setExpectedMomentumAmpDiff->SetGuidance( "\n    Expected difference "
        "between momenta of the beam\n    on its start and in the target. This "
        "parameter can\n    be used to automatically calculate value of the\n"
        "    previous parameter 'momentumAmp'" );
    setExpectedMomentumAmpDiff->SetParameterName( "RecMomentumAmpDiff", false );
    setExpectedMomentumAmpDiff->SetDefaultValue( 0 );
    setExpectedMomentumAmpDiff->SetUnitCandidates( "eV keV MeV GeV" );
    setExpectedMomentumAmpDiff->SetDefaultUnit( "MeV" );
    setExpectedMomentumAmpDiff->AvailableForStates( G4State_PreInit,
                                                    G4State_Idle );

    setEDCollectionAlgorithm = new G4UIcmdWithAString(
        ( CexmcMessenger::reconstructorDirName + "edCollectionAlgo" ).c_str(),
        this );
    setEDCollectionAlgorithm->SetGuidance(
        "\n    Choose crystals to be selected when energy deposit collected\n"
        "    all - all,\n"
        "    adjacent - crystal with maximum energy deposit and\n"
        "               adjacent crystals" );
    setEDCollectionAlgorithm->SetParameterName( "EDCollectionAlgo", false );
    setEDCollectionAlgorithm->SetCandidates( "all adjacent" );
    setEDCollectionAlgorithm->SetDefaultValue( "all" );
    setEDCollectionAlgorithm->AvailableForStates( G4State_PreInit,
                                                  G4State_Idle );
}


CexmcChargeExchangeReconstructorMessenger::
                                    ~CexmcChargeExchangeReconstructorMessenger()
{
    delete useTableMass;
    delete useMassCut;
    delete mCutOPCenter;
    delete mCutNOPCenter;
    delete mCutOPWidth;
    delete mCutNOPWidth;
    delete mCutAngle;
    delete useAbsorbedEnergyCut;
    delete aeCutCLCenter;
    delete aeCutCRCenter;
    delete aeCutCLWidth;
    delete aeCutCRWidth;
    delete aeCutAngle;
    delete setExpectedMomentumAmp;
    delete setExpectedMomentumAmpDiff;
    delete setEDCollectionAlgorithm;
}


void  CexmcChargeExchangeReconstructorMessenger::SetNewValue(
                                        G4UIcommand *  cmd, G4String  value )
{
    do
    {
        if ( cmd == useTableMass )
        {
            reconstructor->UseTableMass(
                        G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
        if ( cmd == useMassCut )
        {
            reconstructor->UseMassCut(
                        G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
        if ( cmd == mCutOPCenter )
        {
            reconstructor->SetMassCutOPCenter(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == mCutNOPCenter )
        {
            reconstructor->SetMassCutNOPCenter(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == mCutOPWidth )
        {
            reconstructor->SetMassCutOPWidth(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == mCutNOPWidth )
        {
            reconstructor->SetMassCutNOPWidth(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == mCutAngle )
        {
            reconstructor->SetMassCutEllipseAngle(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == useAbsorbedEnergyCut )
        {
            reconstructor->UseAbsorbedEnergyCut(
                        G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
        if ( cmd == aeCutCLCenter )
        {
            reconstructor->SetAbsorbedEnergyCutCLCenter(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == aeCutCRCenter )
        {
            reconstructor->SetAbsorbedEnergyCutCRCenter(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == aeCutCLWidth )
        {
            reconstructor->SetAbsorbedEnergyCutCLWidth(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == aeCutCRWidth )
        {
            reconstructor->SetAbsorbedEnergyCutCRWidth(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == aeCutAngle )
        {
            reconstructor->SetAbsorbedEnergyCutEllipseAngle(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setExpectedMomentumAmp )
        {
            reconstructor->SetExpectedMomentumAmp(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setExpectedMomentumAmpDiff )
        {
            reconstructor->SetExpectedMomentumAmpDiff(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setEDCollectionAlgorithm )
        {
            CexmcEDCollectionAlgoritm
                        edCollectionAlgorithm( CexmcCollectEDInAllCrystals );
            do
            {
                if ( value == "adjacent" )
                {
                    edCollectionAlgorithm = CexmcCollectEDInAdjacentCrystals;
                    break;
                }
            } while ( false );
            reconstructor->SetEDCollectionAlgorithm( edCollectionAlgorithm );
            break;
        }
    } while ( false );
}


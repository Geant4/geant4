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
 *       Filename:  CexmcEnergyDepositDigitizerMessenger.cc
 *
 *    Description:  energy deposit digitizer messenger
 *
 *        Version:  1.0
 *        Created:  29.11.2009 19:07:05
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWithoutParameter.hh>
#include "CexmcEnergyDepositDigitizer.hh"
#include "CexmcEnergyDepositDigitizerMessenger.hh"
#include "CexmcMessenger.hh"
#include "CexmcCommon.hh"


CexmcEnergyDepositDigitizerMessenger::CexmcEnergyDepositDigitizerMessenger(
        CexmcEnergyDepositDigitizer *  energyDepositDigitizer_ ) :
    energyDepositDigitizer( energyDepositDigitizer_ ),
    setMonitorThreshold( NULL ), setVetoCountersThreshold( NULL ),
    setLeftVetoCounterThreshold( NULL ), setRightVetoCounterThreshold( NULL ),
    setCalorimetersThreshold( NULL ), setLeftCalorimeterThreshold( NULL ),
    setRightCalorimeterThreshold( NULL ),
    setCalorimeterTriggerAlgorithm( NULL ),
    setOuterCrystalsVetoAlgorithm( NULL ), setOuterCrystalsVetoFraction( NULL ),
    applyFiniteCrystalResolution( NULL ), addCrystalResolutionRange( NULL ),
    clearCrystalResolutionData( NULL )
{
    setMonitorThreshold = new G4UIcmdWithADoubleAndUnit(
            ( CexmcMessenger::monitorEDDirName + "threshold" ).c_str(), this );
    setMonitorThreshold->SetGuidance( "Monitor trigger threshold" );
    setMonitorThreshold->SetParameterName( "MonitorThreshold", false );
    setMonitorThreshold->SetDefaultValue( 0 );
    setMonitorThreshold->SetUnitCandidates( "ev keV MeV GeV" );
    setMonitorThreshold->SetDefaultUnit( "MeV" );
    setMonitorThreshold->AvailableForStates( G4State_PreInit, G4State_Idle );

    setVetoCountersThreshold = new G4UIcmdWithADoubleAndUnit(
            ( CexmcMessenger::vetoCounterEDDirName + "threshold" ).c_str(),
            this );
    setVetoCountersThreshold->SetGuidance( "Veto counters trigger threshold" );
    setVetoCountersThreshold->SetParameterName( "VetoCountersThreshold",
                                                false );
    setVetoCountersThreshold->SetDefaultValue( 0 );
    setVetoCountersThreshold->SetUnitCandidates( "ev keV MeV GeV" );
    setVetoCountersThreshold->SetDefaultUnit( "MeV" );
    setVetoCountersThreshold->AvailableForStates( G4State_PreInit,
                                                  G4State_Idle );

    setLeftVetoCounterThreshold = new G4UIcmdWithADoubleAndUnit(
            ( CexmcMessenger::vetoCounterLeftEDDirName + "threshold" ).c_str(),
            this );
    setLeftVetoCounterThreshold->SetGuidance(
                                    "Left veto counter trigger threshold" );
    setLeftVetoCounterThreshold->SetParameterName( "LeftVetoCounterThreshold",
                                                   false );
    setLeftVetoCounterThreshold->SetDefaultValue( 0 );
    setLeftVetoCounterThreshold->SetUnitCandidates( "ev keV MeV GeV" );
    setLeftVetoCounterThreshold->SetDefaultUnit( "MeV" );
    setLeftVetoCounterThreshold->AvailableForStates( G4State_PreInit,
                                                     G4State_Idle );

    setRightVetoCounterThreshold = new G4UIcmdWithADoubleAndUnit(
            ( CexmcMessenger::vetoCounterRightEDDirName + "threshold" ).c_str(),
            this );
    setRightVetoCounterThreshold->SetGuidance(
                                    "Right veto counter trigger threshold" );
    setRightVetoCounterThreshold->SetParameterName( "RightVetoCounterThreshold",
                                                    false );
    setRightVetoCounterThreshold->SetDefaultValue( 0 );
    setRightVetoCounterThreshold->SetUnitCandidates( "ev keV MeV GeV" );
    setRightVetoCounterThreshold->SetDefaultUnit( "MeV" );
    setRightVetoCounterThreshold->AvailableForStates( G4State_PreInit,
                                                      G4State_Idle );

    setCalorimetersThreshold = new G4UIcmdWithADoubleAndUnit(
            ( CexmcMessenger::calorimeterEDDirName + "threshold" ).c_str(),
            this );
    setCalorimetersThreshold->SetGuidance( "Calorimeters trigger threshold" );
    setCalorimetersThreshold->SetParameterName( "CalorimetersThreshold",
                                                false );
    setCalorimetersThreshold->SetDefaultValue( 0 );
    setCalorimetersThreshold->SetUnitCandidates( "ev keV MeV GeV" );
    setCalorimetersThreshold->SetDefaultUnit( "MeV" );
    setCalorimetersThreshold->AvailableForStates( G4State_PreInit,
                                                  G4State_Idle );

    setLeftCalorimeterThreshold = new G4UIcmdWithADoubleAndUnit(
            ( CexmcMessenger::calorimeterLeftEDDirName + "threshold" ).c_str(),
            this );
    setLeftCalorimeterThreshold->SetGuidance(
                                    "Left calorimeter trigger threshold" );
    setLeftCalorimeterThreshold->SetParameterName( "LeftCalorimeterThreshold",
                                                   false );
    setLeftCalorimeterThreshold->SetDefaultValue( 0 );
    setLeftCalorimeterThreshold->SetUnitCandidates( "ev keV MeV GeV" );
    setLeftCalorimeterThreshold->SetDefaultUnit( "MeV" );
    setLeftCalorimeterThreshold->AvailableForStates( G4State_PreInit,
                                                     G4State_Idle );

    setRightCalorimeterThreshold = new G4UIcmdWithADoubleAndUnit(
            ( CexmcMessenger::calorimeterRightEDDirName + "threshold" ).c_str(),
            this );
    setRightCalorimeterThreshold->SetGuidance(
                                    "Right calorimeter trigger threshold" );
    setRightCalorimeterThreshold->SetParameterName( "RightCalorimeterThreshold",
                                                    false );
    setRightCalorimeterThreshold->SetDefaultValue( 0 );
    setRightCalorimeterThreshold->SetUnitCandidates( "ev keV MeV GeV" );
    setRightCalorimeterThreshold->SetDefaultUnit( "MeV" );
    setRightCalorimeterThreshold->AvailableForStates( G4State_PreInit,
                                                      G4State_Idle );

    setCalorimeterTriggerAlgorithm = new G4UIcmdWithAString(
            ( CexmcMessenger::detectorDirName +
              "calorimeterTriggerAlgorithm" ).c_str(), this );
    setCalorimeterTriggerAlgorithm->SetGuidance( "\n"
                "    all - energy deposit in all crystals in a calorimeter\n"
                "          will be checked against calorimeter threshold "
                          "value,\n"
                "    inner - energy deposit in only inner crystals\n"
                "            will be checked against calorimeter threshold "
                            "value" );
    setCalorimeterTriggerAlgorithm->SetParameterName(
                                        "CalorimeterTriggerAlgorithm", false );
    setCalorimeterTriggerAlgorithm->SetCandidates( "all inner" );
    setCalorimeterTriggerAlgorithm->SetDefaultValue( "inner" );
    setCalorimeterTriggerAlgorithm->AvailableForStates( G4State_PreInit,
                                                        G4State_Idle );

    setOuterCrystalsVetoAlgorithm = new G4UIcmdWithAString(
            ( CexmcMessenger::detectorDirName +
              "outerCrystalsVetoAlgorithm" ).c_str(), this );
    setOuterCrystalsVetoAlgorithm->SetGuidance( "\n"
                "    none - events will not be rejected by any algorithm,\n"
                "    max - reject event trigger if crystal with maximum energy "
                    "\n          deposit is one of outer crystals,\n"
                "    fraction - reject event trigger if energy deposit "
                    "fraction in\n               outer crystals is more than "
                    "value of\n               'outerCrystalsVetoFraction'" );
    setOuterCrystalsVetoAlgorithm->SetParameterName(
                                        "OuterCrystalsVetoAlgorithm", false );
    setOuterCrystalsVetoAlgorithm->SetCandidates( "none max fraction" );
    setOuterCrystalsVetoAlgorithm->SetDefaultValue( "none" );
    setOuterCrystalsVetoAlgorithm->AvailableForStates( G4State_PreInit,
                                                       G4State_Idle );

    setOuterCrystalsVetoFraction = new G4UIcmdWithADouble(
            ( CexmcMessenger::detectorDirName +
              "outerCrystalsVetoFraction" ).c_str(), this );
    setOuterCrystalsVetoFraction->SetGuidance( "\n    Fraction of whole energy "
            "deposit in one calorimeter\n    that belongs to outer crystals.\n"
            "    If 'outerCrystalsVetoAlgorithm' is 'fraction' and\n"
            "    the outer crystals energy deposit fraction exceeds "
                "this\n    value then event won't trigger" );
    setOuterCrystalsVetoFraction->SetParameterName(
                                        "OuterCrystalsVetoFraction", false );
    setOuterCrystalsVetoFraction->SetDefaultValue( 0 );
    setOuterCrystalsVetoFraction->AvailableForStates( G4State_PreInit,
                                                      G4State_Idle );

    applyFiniteCrystalResolution = new G4UIcmdWithABool(
            ( CexmcMessenger::detectorDirName +
              "applyFiniteCrystalResolution" ).c_str(), this );
    applyFiniteCrystalResolution->SetGuidance( "\n     Specify if finite "
            "energy resolution of the crystals\n     will be accounted" );
    applyFiniteCrystalResolution->SetParameterName(
                                        "ApplyFiniteCrystalResolution", true );
    applyFiniteCrystalResolution->SetDefaultValue( true );
    applyFiniteCrystalResolution->AvailableForStates( G4State_PreInit,
                                                      G4State_Idle );

    addCrystalResolutionRange = new G4UIcmdWith3Vector(
            ( CexmcMessenger::detectorDirName +
              "addCrystalResolutionRange" ).c_str(), this );
    addCrystalResolutionRange->SetGuidance( "\n     Add new energy range "
            "(in GeV!) with fwhm percentage\n     value of crystal resolution "
            "in this range" );
    addCrystalResolutionRange->SetParameterName(
            "CrystalResolutionRangeBottom", "CrystalResolutionRangeTop",
            "CrystalResolutionRangeValue", false );
    addCrystalResolutionRange->SetRange( "CrystalResolutionRangeBottom >= 0. "
            "&& CrystalResolutionRangeTop >= 0. && "
            "CrystalResolutionRangeValue >= 0." );
    addCrystalResolutionRange->AvailableForStates( G4State_PreInit,
                                                   G4State_Idle );

    clearCrystalResolutionData = new G4UIcmdWithoutParameter(
            ( CexmcMessenger::detectorDirName +
              "clearCrystalResolutionData" ).c_str(), this );
    clearCrystalResolutionData->SetGuidance( "\n     Clear all crystal "
              "resolution ranges.\n     Can be used to redefine crystal "
              "resolution data" );
    clearCrystalResolutionData->AvailableForStates( G4State_PreInit,
                                                    G4State_Idle );
}


CexmcEnergyDepositDigitizerMessenger::~CexmcEnergyDepositDigitizerMessenger()
{
    delete setMonitorThreshold;
    delete setVetoCountersThreshold;
    delete setLeftVetoCounterThreshold;
    delete setRightVetoCounterThreshold;
    delete setCalorimetersThreshold;
    delete setLeftCalorimeterThreshold;
    delete setRightCalorimeterThreshold;
    delete setCalorimeterTriggerAlgorithm;
    delete setOuterCrystalsVetoAlgorithm;
    delete setOuterCrystalsVetoFraction;
    delete applyFiniteCrystalResolution;
    delete addCrystalResolutionRange;
    delete clearCrystalResolutionData;
}


void  CexmcEnergyDepositDigitizerMessenger::SetNewValue( G4UIcommand *  cmd,
                                                         G4String  value )
{
    do
    {
        if ( cmd == setMonitorThreshold )
        {
            energyDepositDigitizer->SetMonitorThreshold(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setVetoCountersThreshold )
        {
            energyDepositDigitizer->SetVetoCountersThreshold(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setLeftVetoCounterThreshold )
        {
            energyDepositDigitizer->SetVetoCounterLeftThreshold(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setRightVetoCounterThreshold )
        {
            energyDepositDigitizer->SetVetoCounterRightThreshold(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setCalorimetersThreshold )
        {
            energyDepositDigitizer->SetCalorimetersThreshold(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setLeftCalorimeterThreshold )
        {
            energyDepositDigitizer->SetCalorimeterLeftThreshold(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setRightCalorimeterThreshold )
        {
            energyDepositDigitizer->SetCalorimeterRightThreshold(
                        G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == setCalorimeterTriggerAlgorithm )
        {
            CexmcCalorimeterTriggerAlgorithm  calorimeterTriggerAlgorithm(
                                    CexmcAllCrystalsMakeEDTriggerThreshold );
            do
            {
                if ( value == "inner" )
                {
                    calorimeterTriggerAlgorithm =
                                    CexmcInnerCrystalsMakeEDTriggerThreshold;
                    break;
                }
            } while ( false );
            energyDepositDigitizer->SetCalorimeterTriggerAlgorithm(
                                                calorimeterTriggerAlgorithm );
            break;
        }
        if ( cmd == setOuterCrystalsVetoAlgorithm )
        {
            CexmcOuterCrystalsVetoAlgorithm  outerCrystalsVetoAlgorithm(
                                                    CexmcNoOuterCrystalsVeto );
            do
            {
                if ( value == "max" )
                {
                    outerCrystalsVetoAlgorithm =
                                        CexmcMaximumEDInASingleOuterCrystalVeto;
                    break;
                }
                if ( value == "fraction" )
                {
                    outerCrystalsVetoAlgorithm =
                                        CexmcFractionOfEDInOuterCrystalsVeto;
                    break;
                }
            } while ( false );
            energyDepositDigitizer->SetOuterCrystalsVetoAlgorithm(
                                                outerCrystalsVetoAlgorithm );
            break;
        }
        if ( cmd == setOuterCrystalsVetoFraction )
        {
            energyDepositDigitizer->SetOuterCrystalsVetoFraction(
                        G4UIcmdWithADouble::GetNewDoubleValue( value ) );
            break;
        }
        if ( cmd == applyFiniteCrystalResolution )
        {
            energyDepositDigitizer->ApplyFiniteCrystalResolution(
                        G4UIcmdWithABool::GetNewBoolValue( value ) );
            break;
        }
        if ( cmd == addCrystalResolutionRange )
        {
            G4ThreeVector  vec( G4UIcmdWith3Vector::GetNew3VectorValue(
                                                                    value ) );
            G4double       bottom( std::min( vec.x(), vec.y() ) );
            G4double       top( std::max( vec.x(), vec.y() ) );
            energyDepositDigitizer->AddCrystalResolutionRange( bottom, top,
                                                               vec.z() );
            break;
        }
        if ( cmd == clearCrystalResolutionData )
        {
            energyDepositDigitizer->ClearCrystalResolutionData();
            break;
        }
    } while ( false );
}


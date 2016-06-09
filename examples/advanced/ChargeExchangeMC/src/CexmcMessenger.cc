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
 *       Filename:  CexmcMessenger.cc
 *
 *    Description:  common messenger stuff (directories etc.)
 *
 *        Version:  1.0
 *        Created:  15.11.2009 12:59:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4UIdirectory.hh>
#include "CexmcMessenger.hh"
#include "CexmcSensitiveDetectorsAttributes.hh"


namespace
{
    const G4String  ed( CexmcDetectorTypeName[ CexmcEDDetector ] + "/" );
}


CexmcMessenger *  CexmcMessenger::instance( NULL );

G4String  CexmcMessenger::mainDirName( "/cexmc/" );
G4String  CexmcMessenger::geometryDirName( CexmcMessenger::mainDirName +
                                           "geometry/" );
G4String  CexmcMessenger::physicsDirName( CexmcMessenger::mainDirName +
                                           "physics/" );
G4String  CexmcMessenger::gunDirName( CexmcMessenger::mainDirName +
                                           "gun/" );
G4String  CexmcMessenger::detectorDirName( CexmcMessenger::mainDirName +
                                           "detector/" );
G4String  CexmcMessenger::eventDirName( CexmcMessenger::mainDirName +
                                           "event/" );
G4String  CexmcMessenger::runDirName( CexmcMessenger::mainDirName +
                                           "run/" );
G4String  CexmcMessenger::monitorDirName( CexmcMessenger::detectorDirName +
               CexmcDetectorRoleName[ CexmcMonitorDetectorRole ] + "/" );
G4String  CexmcMessenger::targetDirName( CexmcMessenger::detectorDirName +
               CexmcDetectorRoleName[ CexmcTargetDetectorRole ] + "/" );
G4String  CexmcMessenger::vetoCounterDirName( CexmcMessenger::detectorDirName +
               CexmcDetectorRoleName[ CexmcVetoCounterDetectorRole ] + "/" );
G4String  CexmcMessenger::vetoCounterLeftDirName(
                            CexmcMessenger::vetoCounterDirName + "left/" );
G4String  CexmcMessenger::vetoCounterRightDirName(
                            CexmcMessenger::vetoCounterDirName + "right/" );
G4String  CexmcMessenger::calorimeterDirName( CexmcMessenger::detectorDirName +
               CexmcDetectorRoleName[ CexmcCalorimeterDetectorRole ] + "/" );
G4String  CexmcMessenger::calorimeterLeftDirName(
                            CexmcMessenger::calorimeterDirName + "left/" );
G4String  CexmcMessenger::calorimeterRightDirName(
                            CexmcMessenger::calorimeterDirName + "right/" );
G4String  CexmcMessenger::monitorEDDirName(
                            CexmcMessenger::monitorDirName + ed );
G4String  CexmcMessenger::vetoCounterEDDirName(
                            CexmcMessenger::vetoCounterDirName + ed );
G4String  CexmcMessenger::vetoCounterLeftEDDirName(
                            CexmcMessenger::vetoCounterLeftDirName + ed );
G4String  CexmcMessenger::vetoCounterRightEDDirName(
                            CexmcMessenger::vetoCounterRightDirName + ed );
G4String  CexmcMessenger::calorimeterEDDirName(
                            CexmcMessenger::calorimeterDirName + ed );
G4String  CexmcMessenger::calorimeterLeftEDDirName(
                            CexmcMessenger::calorimeterLeftDirName + ed );
G4String  CexmcMessenger::calorimeterRightEDDirName(
                            CexmcMessenger::calorimeterRightDirName + ed );
G4String  CexmcMessenger::reconstructorDirName(
                            CexmcMessenger::mainDirName + "reconstructor/" );
G4String  CexmcMessenger::visDirName( CexmcMessenger::mainDirName + "vis/" );
#ifdef CEXMC_USE_ROOT
G4String  CexmcMessenger::histoDirName(
                            CexmcMessenger::mainDirName + "histo/" );
#endif


CexmcMessenger *  CexmcMessenger::Instance( void )
{
    if ( instance == NULL )
        instance = new CexmcMessenger;

    return instance;
}


void  CexmcMessenger::Destroy( void )
{
    delete instance;
    instance = NULL;
}


CexmcMessenger::CexmcMessenger() : mainDir( NULL ), geometryDir( NULL ),
    physicsDir( NULL ), gunDir( NULL ), detectorDir( NULL ), eventDir( NULL ),
    runDir( NULL ), monitorDir( NULL ), targetDir( NULL ),
    vetoCounterDir( NULL ), vetoCounterLeftDir( NULL ),
    vetoCounterRightDir( NULL ), calorimeterDir( NULL ),
    calorimeterLeftDir( NULL ), calorimeterRightDir( NULL ),
    monitorEDDir( NULL ), vetoCounterEDDir( NULL ),
    vetoCounterLeftEDDir( NULL ), vetoCounterRightEDDir( NULL ),
    calorimeterEDDir( NULL ), calorimeterLeftEDDir( NULL ),
    calorimeterRightEDDir( NULL ), reconstructorDir( NULL ), visDir( NULL )
#ifdef CEXMC_USE_ROOT
    ,histoDir( NULL )
#endif
{
    mainDir = new G4UIdirectory( mainDirName );
    mainDir->SetGuidance( "Cexmc settings." );
    geometryDir = new G4UIdirectory( geometryDirName );
    geometryDir->SetGuidance( "Geometry settings (gdml file etc.)" );
    physicsDir = new G4UIdirectory( physicsDirName );
    physicsDir->SetGuidance( "Physics related settings "
                             "(production model etc.)" );
    gunDir = new G4UIdirectory( gunDirName );
    gunDir->SetGuidance( "Gun settings (different FWHMs etc.)" );
    detectorDir = new G4UIdirectory( detectorDirName );
    detectorDir->SetGuidance( "Sensitive detectors settings" );
    eventDir = new G4UIdirectory( eventDirName );
    eventDir->SetGuidance( "Event settings (verbose level etc.)" );
    runDir = new G4UIdirectory( runDirName );
    runDir->SetGuidance( "Run settings (geometry file etc.)" );
    monitorDir = new G4UIdirectory( monitorDirName );
    monitorDir->SetGuidance( "Various settings for the monitor." );
    targetDir = new G4UIdirectory( targetDirName );
    targetDir->SetGuidance( "Various settings for the target." );
    vetoCounterDir = new G4UIdirectory( vetoCounterDirName );
    vetoCounterDir->SetGuidance( "Various settings for the veto counters." );
    vetoCounterLeftDir = new G4UIdirectory( vetoCounterLeftDirName );
    vetoCounterLeftDir->SetGuidance(
            "\n    Various settings for the left veto counter." );
    vetoCounterRightDir = new G4UIdirectory( vetoCounterRightDirName );
    vetoCounterRightDir->SetGuidance(
            "\n    Various settings for the right veto counter." );
    calorimeterDir = new G4UIdirectory( calorimeterDirName );
    calorimeterDir->SetGuidance( "Various settings for the calorimeters." );
    calorimeterLeftDir = new G4UIdirectory( calorimeterLeftDirName );
    calorimeterLeftDir->SetGuidance(
            "\n    Various settings for the left calorimeter." );
    calorimeterRightDir = new G4UIdirectory( calorimeterRightDirName );
    calorimeterRightDir->SetGuidance(
            "\n    Various settings for the right calorimeter." );
    monitorEDDir = new G4UIdirectory( monitorEDDirName );
    monitorEDDir->SetGuidance(
            "\n    Energy deposit settings in the monitor (thresholds etc.)" );
    vetoCounterEDDir = new G4UIdirectory( vetoCounterEDDirName );
    vetoCounterEDDir->SetGuidance(
            "\n    Energy deposit settings in both veto counters "
            "(thresholds etc.)" );
    vetoCounterLeftEDDir = new G4UIdirectory( vetoCounterLeftEDDirName );
    vetoCounterLeftEDDir->SetGuidance(
            "\n    Energy deposit settings in the left veto counter "
            "(thresholds etc.)" );
    vetoCounterRightEDDir = new G4UIdirectory( vetoCounterRightEDDirName );
    vetoCounterRightEDDir->SetGuidance(
            "\n    Energy deposit settings in the right veto counter "
            "(thresholds etc.)" );
    calorimeterEDDir = new G4UIdirectory( calorimeterEDDirName );
    calorimeterEDDir->SetGuidance(
            "\n    Energy deposit settings in both calorimeters "
            "(thresholds etc.)" );
    calorimeterLeftEDDir = new G4UIdirectory( calorimeterLeftEDDirName );
    calorimeterLeftEDDir->SetGuidance(
            "\n    Energy deposit settings in the left calorimeter "
            "(thresholds etc.)" );
    calorimeterRightEDDir = new G4UIdirectory( calorimeterRightEDDirName );
    calorimeterRightEDDir->SetGuidance(
            "\n    Energy deposit settings in the right calorimeter "
            "(thresholds etc.)" );
    reconstructorDir = new G4UIdirectory( reconstructorDirName );
    reconstructorDir->SetGuidance( "Reconstructor settings" );
    visDir = new G4UIdirectory( visDirName );
    visDir->SetGuidance( "Visualization settings" );
#ifdef CEXMC_USE_ROOT
    histoDir = new G4UIdirectory( histoDirName );
    histoDir->SetGuidance( "Commands to list and show histograms" );
#endif
}


CexmcMessenger::~CexmcMessenger()
{
    delete mainDir;
    delete geometryDir;
    delete physicsDir;
    delete gunDir;
    delete detectorDir;
    delete eventDir;
    delete runDir;
    delete monitorDir;
    delete targetDir;
    delete vetoCounterDir;
    delete vetoCounterLeftDir;
    delete vetoCounterRightDir;
    delete calorimeterDir;
    delete calorimeterLeftDir;
    delete calorimeterRightDir;
    delete monitorEDDir;
    delete vetoCounterEDDir;
    delete vetoCounterLeftEDDir;
    delete vetoCounterRightEDDir;
    delete calorimeterEDDir;
    delete calorimeterLeftEDDir;
    delete calorimeterRightEDDir;
    delete reconstructorDir;
    delete visDir;
#ifdef CEXMC_USE_ROOT
    delete histoDir;
#endif
}


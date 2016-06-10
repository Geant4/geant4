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
//
// $Id: G4ApplicationState.hh 67970 2013-03-13 10:10:06Z gcosmo $
//

#ifndef G4APPLICATIONSTATE_H
#define G4APPLICATIONSTATE_H 1

// Class Description:
//
// Specifies the state of the G4 application
//
// States:
//  G4State_PreInit
//    At the very begining of the Application. G4StateManager starts
//    with this state. G4Initializer changes this state to Init when
//    G4Initializer::Initialize() method starts. At the moment of
//    the state change of PreInit->Init, no material, geometrical,
//    particle or physics process has been initialized.
//  G4State_Init
//    During the G4Initializer::Initialize() method. G4Initializer
//    changes this state to Idle when all initialization procedures
//    are successfully Done.
//  G4State_Idle
//    Ready to start "BeamOn". G4RunManager changes this state to
//    GeomClosed when G4RunManager::BeamOn() method starts and
//    G4GeometryManager::CloseGeometry() is Done. At the end of
//    BeamOn() method, G4RunManager will reset the application state
//    to Idle after G4GeometryManager::OpenGeometry() is Done.
//  G4State_GeomClosed
//    G4 is in this state between G4GeometryManager::CloseGeometry()
//    and G4GeometryManager::OpenGeometry(), but no event is in 
//    progress. At the begining of each event (construction of a
//    G4Event object and primary particle generation), G4RunManager
//    changes this state to EventProc and resets to GeomClosed state
//    when G4EventManager::ProcessOneEvent() is over.
//  G4State_EventProc
//    Processing an event.
//  G4State_Quit
//    G4 is in this state when the destructor of G4RunManager is invoked.
//  G4State_Abort
//    G4 is in this state when G4Exception is invoked.
//    
//
//  PreInit
//    |
//    v
//  Init
//    |
//    v
//  Idle ------> Quit
//    |^
//    v|
//  GeomClosed (at each run)
//    |^
//    v|
//  EventProc (at each event)
//

enum G4ApplicationState 
  {G4State_PreInit, G4State_Init, G4State_Idle, G4State_GeomClosed,
   G4State_EventProc, G4State_Quit, G4State_Abort};

#endif


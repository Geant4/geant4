//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ApplicationState.hh,v 1.2 2001-07-11 10:00:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4APPLICATIONSTATE_H
#define G4APPLICATIONSTATE_H 1

// Class Description:
//
// Specifies the state of the G4 application
//
// States:
//  PreInit
//    At the very begining of the Application. G4StateManager starts
//    with this state. G4Initializer changes this state to Init when
//    G4Initializer::Initialize() method starts. At the moment of
//    the state change of PreInit->Init, no material, geometrical,
//    particle or physics process has been initialized.
//  Init
//    During the G4Initializer::Initialize() method. G4Initializer
//    changes this state to Idle when all initialization procedures
//    are successfully Done.
//  Idle
//    Ready to start "BeamOn". G4RunManager changes this state to
//    GeomClosed when G4RunManager::BeamOn() method starts and
//    G4GeometryManager::CloseGeometry() is Done. At the end of
//    BeamOn() method, G4RunManager will reset the application state
//    to Idle after G4GeometryManager::OpenGeometry() is Done.
//  GeomClosed
//    G4 is in this state between G4GeometryManager::CloseGeometry()
//    and G4GeometryManager::OpenGeometry(), but no event is in 
//    progress. At the begining of each event (construction of a
//    G4Event object and primary particle generation), G4RunManager
//    changes this state to EventProc and resets to GeomClosed state
//    when G4EventManager::ProcessOneEvent() is over.
//  EventProc
//    Processing an event.
//  Quit
//    G4 is in this state when the destructor of G4RunManager is invoked.
//  Abort
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
  {PreInit, Init, Idle, GeomClosed, EventProc, Quit, Abort};

#endif


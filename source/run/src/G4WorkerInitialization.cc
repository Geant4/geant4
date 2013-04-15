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
// $Id$
//
// Author: Ivana Hrivnacova, 10/04/2013  (ivana@ipno.in2p3.fr)

#include "G4WorkerInitialization.hh"
#include "G4VUserApplication.hh"

//_____________________________________________________________________________
G4WorkerInitialization::G4WorkerInitialization(G4VUserApplication* userApplication)
 : G4VUserWorkerInitialization(),
   fUserApplication(userApplication)
{}

//_____________________________________________________________________________
G4WorkerInitialization::~G4WorkerInitialization()
{}  

//_____________________________________________________________________________
void G4WorkerInitialization::WorkerStart() const
{
  //G4cout << "G4WorkerInitialization::WorkerStart " << this << G4endl;

  // Create new user application instance for the worker
  // (so that it handles only the user classes created on the worker) 
  G4VUserApplication* userApplication = fUserApplication->CreateInstance();

  SetUserAction(userApplication->CreatePrimaryGeneratorAction());
  
  G4UserRunAction* runAction = userApplication->CreateRunAction();
  //G4cout << "runAction " << runAction << G4endl; 
  if ( runAction) SetUserAction(runAction);
  
  G4UserEventAction* eventAction = userApplication->CreateEventAction();
  //G4cout << "eventAction " << eventAction << G4endl; 
  if (eventAction) SetUserAction(eventAction);
  
  G4UserTrackingAction* trackingAction = userApplication->CreateTrackingAction();
  //G4cout << "trackingAction " << trackingAction << G4endl; 
  if (trackingAction) SetUserAction(trackingAction);
  
  G4UserSteppingAction* steppingAction = userApplication->CreateSteppingAction();
  //G4cout << "steppingAction " << steppingAction << G4endl; 
  if (steppingAction) SetUserAction(steppingAction);
  
  G4UserStackingAction* stackingAction = userApplication->CreateStackingAction();
  //G4cout << "stackingAction " << stackingAction << G4endl; 
  if (stackingAction) SetUserAction(stackingAction);
  
  // Perform additional setting (if needed)
  userApplication->Initialize();
    
  // Delete the worker user application
  delete userApplication;

  //G4cout << "G4WorkerInitialization::WorkerStart done " << this << G4endl;
}  


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
// $Id: $
//

#include "G4RTWorkerInitialization.hh"
#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserStackingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4UserSteppingAction.hh"

#include "G4RTRunAction.hh"
#include "G4RTPrimaryGeneratorAction.hh"
#include "G4RTTrackingAction.hh"
#include "G4RTSteppingAction.hh"

G4ThreadLocal const G4UserRunAction * G4RTWorkerInitialization::theUserRunAction = 0;
G4ThreadLocal const G4VUserPrimaryGeneratorAction * G4RTWorkerInitialization::theUserPrimaryGeneratorAction = 0;
G4ThreadLocal const G4UserEventAction * G4RTWorkerInitialization::theUserEventAction = 0;
G4ThreadLocal const G4UserStackingAction * G4RTWorkerInitialization::theUserStackingAction = 0;
G4ThreadLocal const G4UserTrackingAction * G4RTWorkerInitialization::theUserTrackingAction = 0;
G4ThreadLocal const G4UserSteppingAction * G4RTWorkerInitialization::theUserSteppingAction = 0;

G4ThreadLocal G4RTRunAction * G4RTWorkerInitialization::theRTRunAction = 0;
G4ThreadLocal G4RTPrimaryGeneratorAction * G4RTWorkerInitialization::theRTPrimaryGeneratorAction = 0;
G4ThreadLocal G4RTTrackingAction * G4RTWorkerInitialization::theRTTrackingAction = 0;
G4ThreadLocal G4RTSteppingAction * G4RTWorkerInitialization::theRTSteppingAction = 0;

#include "G4WorkerRunManager.hh"

G4RTWorkerInitialization::G4RTWorkerInitialization()
{;}

G4RTWorkerInitialization::~G4RTWorkerInitialization()
{;}

void G4RTWorkerInitialization::WorkerRunStart() const
{
  if(!theRTRunAction) theRTRunAction = new G4RTRunAction;
  if(!theRTPrimaryGeneratorAction) theRTPrimaryGeneratorAction = new G4RTPrimaryGeneratorAction;
  if(!theRTTrackingAction) theRTTrackingAction = new G4RTTrackingAction;
  if(!theRTSteppingAction) theRTSteppingAction = new G4RTSteppingAction;
  
  G4WorkerRunManager* runMan = G4WorkerRunManager::GetWorkerRunManager(); 

  theUserRunAction = runMan->GetUserRunAction();
  theUserPrimaryGeneratorAction = runMan->GetUserPrimaryGeneratorAction();
  theUserEventAction = runMan->GetUserEventAction();
  theUserStackingAction = runMan->GetUserStackingAction();
  theUserTrackingAction = runMan->GetUserTrackingAction();
  theUserSteppingAction = runMan->GetUserSteppingAction();
  
  runMan->SetUserAction(theRTRunAction);  
  runMan->SetUserAction(theRTPrimaryGeneratorAction);  
  runMan->SetUserAction(static_cast<G4UserEventAction*>(0));  
  runMan->SetUserAction(static_cast<G4UserStackingAction*>(0));  
  runMan->SetUserAction(theRTTrackingAction);  
  runMan->SetUserAction(theRTSteppingAction);  

  theRTPrimaryGeneratorAction->SetUp();
}

void G4RTWorkerInitialization::WorkerRunEnd() const
{
  G4WorkerRunManager* runMan = G4WorkerRunManager::GetWorkerRunManager(); 
  runMan->SetUserAction(const_cast<G4UserRunAction*>(theUserRunAction));  
  runMan->SetUserAction(const_cast<G4VUserPrimaryGeneratorAction*>(theUserPrimaryGeneratorAction));  
  runMan->SetUserAction(const_cast<G4UserEventAction*>(theUserEventAction));  
  runMan->SetUserAction(const_cast<G4UserStackingAction*>(theUserStackingAction));  
  runMan->SetUserAction(const_cast<G4UserTrackingAction*>(theUserTrackingAction));  
  runMan->SetUserAction(const_cast<G4UserSteppingAction*>(theUserSteppingAction));  
}


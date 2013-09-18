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
#include "G4VUserActionInitialization.hh"
#include "G4RunManager.hh"

G4VUserActionInitialization::G4VUserActionInitialization()
{;}

G4VUserActionInitialization::~G4VUserActionInitialization()
{;}

void G4VUserActionInitialization::BuildForMaster() const
{;}

void G4VUserActionInitialization::SetUserAction(G4VUserPrimaryGeneratorAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4VUserActionInitialization::SetUserAction(G4UserRunAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4VUserActionInitialization::SetUserAction(G4UserEventAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4VUserActionInitialization::SetUserAction(G4UserStackingAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4VUserActionInitialization::SetUserAction(G4UserTrackingAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4VUserActionInitialization::SetUserAction(G4UserSteppingAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

G4VSteppingVerbose* G4VUserActionInitialization::InitializeSteppingVerbose() const
{ return static_cast<G4VSteppingVerbose*>(0); }



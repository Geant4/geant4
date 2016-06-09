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
// $Id: G3toG4RunAction.cc,v 1.4 2006-06-29 17:20:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ios.hh"
#include "G3toG4RunAction.hh"
#include "G4Run.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

G3toG4RunAction::G3toG4RunAction(){
  runIDcounter = 0;
}

G3toG4RunAction::~G3toG4RunAction(){;}

void G3toG4RunAction::BeginOfRunAction(const G4Run* aRun){
  ((G4Run *)(aRun))->SetRunID(runIDcounter++);
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  G4UImanager* UI = G4UImanager::GetUIpointer(); 
  
  if (G4VVisManager::GetConcreteInstance()) {
    UI->ApplyCommand("/vis~/clear/view");
    UI->ApplyCommand("/vis~/draw/current");
  } 
}

void G3toG4RunAction::EndOfRunAction(const G4Run*){;}


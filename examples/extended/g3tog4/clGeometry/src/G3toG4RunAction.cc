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
// $Id: G3toG4RunAction.cc,v 1.2 2001-07-11 09:58:10 gunter Exp $
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

void G3toG4RunAction::EndOfRunAction(const G4Run* aRun){;}


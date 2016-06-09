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
// $Id: perspective.cc,v 1.3 2010-11-09 09:49:57 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "globals.hh"
#include "G4VisExecutive.hh"
#include "G4VisExtent.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

#include "PerspectiveVisAction.hh"

int main(int argc,char** argv) {

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialise();

  visManager->SetUserAction
    (new PerspectiveVisAction,
     G4VisExtent(-5*m,5*m,-5*m,5*m,-5*m,5*m));  // 2nd argument optional.

  G4String command = "/control/execute ";
  G4String fileName = "vis.mac";
  if (argc > 1) {  // Macro name supplied.
    fileName = argv[1];
  }
  G4UIExecutive* ui = new G4UIExecutive(argc,argv);
  G4UImanager::GetUIpointer()->ApplyCommand(command+fileName);    
  ui->SessionStart();

  delete ui;
  delete visManager;
}

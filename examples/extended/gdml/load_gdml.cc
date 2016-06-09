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
// $Id: load_gdml.cc,v 1.6 2007/12/11 13:13:25 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//
// --------------------------------------------------------------
//      GEANT 4 - load_gdml
//
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4TransportationManager.hh"
#include "G4VisExecutive.hh"

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4GDMLParser.hh"

int main(int argc, char **argv) {

   if (argc<2) {
   
      G4cout << "Usage: load_gdml <filename>" << G4endl;
      return 0;
   }

   G4GDMLParser parser;

   parser.Read(argv[1]);

   G4RunManager* runManager = new G4RunManager;
   G4VisManager* visManager = new G4VisExecutive;

   runManager->SetUserInitialization(new DetectorConstruction(parser.GetWorldVolume()));
   runManager->SetUserInitialization(new PhysicsList);
   runManager->SetUserAction(new PrimaryGeneratorAction);

   runManager->Initialize();
   visManager->Initialize();

   G4UImanager* UI = G4UImanager::GetUIpointer();

   G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
   session = new G4UIterminal(new G4UItcsh);
#else
   session = new G4UIterminal();
#endif

   UI->ApplyCommand("/control/execute vis.mac"); 

   session->SessionStart();
   
   delete visManager;
   delete runManager;

   return 0;
}

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
// $Id: hepmcEx02.cc,v 1.4 2006/06/29 17:09:29 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - example of HepMC-interface
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "H02DetectorConstruction.hh"
#include "H02PhysicsList.hh"
#include "H02PrimaryGeneratorAction.hh"
#include "H02EventAction.hh"
#include "H02SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc, char** argv)
{
  G4RunManager* runManager= new G4RunManager;

  runManager-> SetUserInitialization(new H02DetectorConstruction);
  runManager-> SetUserInitialization(new H02PhysicsList);

  runManager-> Initialize();

  runManager-> SetUserAction(new H02PrimaryGeneratorAction);
  runManager-> SetUserAction(new H02EventAction);
  runManager-> SetUserAction(new H02SteppingAction);

#ifdef G4VIS_USE
  // initialize visualization package
  G4VisManager* visManager= new G4VisExecutive;
  visManager-> Initialize();
  G4cout << G4endl;
#endif

  if(argc==1) {
    // G4UIterminal is a (dumb) terminal.

#ifdef QERAUOY
    G4UItcsh* tcsh= new
      G4UItcsh("[40;01;36mhepmcEx02[40;33m(%s)[40;32m[%/][00;30m:");
    G4UIterminal* session= new G4UIterminal(tcsh);
    tcsh-> SetLsColor(RED, GREEN);
#else
    G4UItcsh* tcsh= new G4UItcsh("hepmcEx01(%s)[%/]:");
    G4UIterminal* session= new G4UIterminal(tcsh);
#endif
    session->SessionStart();
    delete session;

  } else {
    G4UImanager* UImanager= G4UImanager::GetUIpointer();
    G4String command= "/control/execute ";
    G4String fileName= argv[1];
    UImanager-> ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}


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
//

#include "Tst12DetectorConstruction.hh"
#include "Tst12RunAction.hh"
#include "Tst12PrimaryGeneratorAction.hh"
#include "Tst12SteppingAction.hh"
#include "Tst12StackingAction.hh"

#include "FTFP_BERT.hh"
#include "FTF_BIC.hh"
#include "QGS_BIC.hh"
#include "QGSP_BERT.hh"
#include "QGSP_FTFP_BERT.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  CLHEP::RanecuEngine defaultEngine;
  CLHEP::HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;
  G4VUserPhysicsList * thePL(0);
  G4String inputFileName="-";
  if (argc > 2) { // second arg is PhysicsList
     G4String opt = argv[2];
     if (opt == "FTFP_BERT" ) thePL=new FTFP_BERT;
     else if (opt == "FTF_BIC" ) thePL=new FTF_BIC;
     else if (opt == "QGS_BIC" ) thePL=new QGS_BIC;
     else if (opt == "QGSP_BERT" ) thePL=new QGSP_BERT;
     else if (opt == "QGSP_FTFP_BERT" ) thePL=new QGSP_FTFP_BERT;

     if ( thePL ) G4cout << "Test12 using physics list " << opt << G4endl;
  }
  if ( !thePL ) {
     thePL=new FTFP_BERT;
     G4cout << "Test12 using physics list FTFP_BERT (default)" << G4endl;
  }
  
  if (argc>1 ) inputFileName=argv[1];
     
  
  // UserInitialization classes
  runManager->SetUserInitialization(new Tst12DetectorConstruction);
  //  runManager->SetUserInitialization(new Tst12PhysicsList);
  runManager->SetUserInitialization(thePL);

  // UserAction classes
  runManager->SetUserAction(new Tst12RunAction);
  runManager->SetUserAction(new Tst12PrimaryGeneratorAction);
  runManager->SetUserAction(new Tst12StackingAction);

  if(inputFileName == "-")
  {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+inputFileName);
  }

  delete runManager;
  return 0;
}


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
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31 ----------------------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"
#include "G4UItcsh.hh"

#include "test31DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "test31PrimaryGeneratorAction.hh"
#include "test31EventAction.hh"
#include "test31TrackingAction.hh"
#include "test31RunAction.hh"
#include "StackingAction.hh"

#include "G4Timer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv) {

  G4Timer* timer = new G4Timer();
  timer->Start();

  G4int verbose = 1;
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  // set mandatory initialization classes
  test31DetectorConstruction* det = new test31DetectorConstruction();
  runManager->SetUserInitialization(det);
  det->SetVerbose(verbose);
  if(verbose >0) G4cout << "Detector Construction is defined" << G4endl;
  
  runManager->SetUserInitialization(new PhysicsList());
  if(verbose >0) G4cout << "Physics List is defined" << G4endl;

  G4cout << "User actions will be initialized" << G4endl;
 
  // set user action classes
  runManager->SetUserAction(new test31PrimaryGeneratorAction());
  if(verbose >0) G4cout << "test31PrimaryGeneratorAction is defined" << G4endl;

  runManager->SetUserAction(new test31RunAction());
  if(verbose >0) G4cout << "test31RunAction is defined" << G4endl;

  test31EventAction* event = new test31EventAction(det);
  runManager->SetUserAction(event);
  if(verbose >0) G4cout << "EventAction is defined" << G4endl;

  runManager->SetUserAction(new StackingAction());
  runManager->SetUserAction(new test31TrackingAction());
  if(verbose >0) G4cout << "TrackingAction is defined" << G4endl;

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(1 < verbose) UI->ListCommands("/test31/");
     
  if (argc==1)   // Define UI terminal for interactive mode
    {
      G4UIsession * session;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);
#else
      session = new G4UIterminal();
#endif
      session->SessionStart();
      delete session;
    }
  else if (argc>1) // Batch mode with 1 or more files
    {
      if(verbose >0) G4cout << "UI interface is started" << G4endl;
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }
    
 // job termination
  
  timer->Stop();
  G4cout << "  "  << *timer << G4endl;
  /*
  G4double R1 = 55.*mm;
  G4double R2 = 400.*mm;
  G4double r  = 82.*mm;
  G4double dr = 15.5*mm;
  G4double lr = std::log(R2/R1);
  G4double al = 1.0;
  G4double d  = R2 - R1;
  G4double dd = d;
  if(al != 1.0) dd = std::pow(R2,2.0 - al) - std::pow(R2,2.0 - al);

  G4cout << "  R1(mm)= " << R1 << "  R2(mm)= " << R2 << G4endl; 
  G4double x;
  for(G4int i=1; i<21; i++) {
    if(al == 1.0)
      x = 2.0*d*(1.0/dd - 1.0/(r*lr));
    else
      x = 2.0*d*(std::pow(r,1.0 - al)/dd - 1.0/((2.0 - al)*r*lr));

    G4cout << "Row# " << i << " r(mm)= " << r << "   E= " << x << G4endl;
    r += dr;
  }
  */
  delete timer;
  delete runManager;
  return 0;
}


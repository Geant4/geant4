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
// $Id: exampleP02.cc,v 1.1 2007/11/29 17:05:21 witoldp Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExP02DetectorConstruction.hh"
#include "ExP02DetConstrReader.hh"
#include "ExP02PhysicsList.hh"
#include "ExP02PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{  

  if(argc==1)
  {
    std::cout << "Please give 'write' or 'read' as argument " << std::endl;
    return 0;
  }

  G4VUserDetectorConstruction* det;

  if(std::string(argv[1]) == "read")
    {
      det = new ExP02DetConstrReader;
    }
  else if(std::string(argv[1]) == "write")
  {
    det = new ExP02DetectorConstruction;
  }
  else
  {
    std::cout << "Wrong argument!" << std::endl;
    return 0;
  }

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(new ExP02PhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new ExP02PrimaryGeneratorAction());

  //Initialize G4 kernel
  runManager->Initialize();
      
  //get the pointer to the User Interface manager 
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  // Define (G)UI terminal for interactive mode  
  // G4UIterminal is a (dumb) terminal.
  G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
  session = new G4UIterminal(new G4UItcsh);      
#else
  session = new G4UIterminal();
#endif    

  UI->ApplyCommand("/control/execute vis.mac");    
  session->SessionStart();
  delete session;
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


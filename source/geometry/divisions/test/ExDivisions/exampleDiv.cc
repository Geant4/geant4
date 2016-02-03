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
// $Id$
//

#include "ExDivDetectorConstruction.hh"
#include "ExDivActionInitialization.hh"

#include "ExDivPhysicsList.hh"
#include "ExDivSteppingVerbose.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "ExDivTesterBox.hh"
#include "ExDivTesterTubs.hh"
#include "ExDivTesterCons.hh"
#include "ExDivTesterTrd.hh"
#include "ExDivTesterPara.hh"
#include "ExDivTesterPolycone.hh"
#include "ExDivTesterPolyhedra.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  G4String theSolidTypeStr="box";
  G4String thePVTypeStr="division";
  G4String thePosTypeStr="normal";
  std::vector<G4String> theExtraPars;

  // First argument is the type of divisioning replica or division
  // the second argument is the type of solid
  std::vector<std::string> vsarg;
  for( G4int jj = 0; jj < argc; jj ++ )
  {
    vsarg.push_back( std::string(argv[jj] ) );
  } 
  
  G4int narg = vsarg.size();
  if( narg == 1 )
  {
    G4cout << "!!! No input division type provided. Defaulting to 'division' "
           << G4endl;
    G4cout << "!!! No positioning type provided. Defaulting to 'normal' "
           << G4endl;
    G4cout << "!!! No input solid type provided. Defaulting to 'box' "
           << G4endl;
  }
  else
  {
    if( narg == 2 )
    {
      G4cout << "!!! No positioning type provided. Defaulting to 'normal' "
             << G4endl;
      G4cout << "!!! No input solid type provided. Defaulting to 'box' "
             << G4endl;
    }
    else
    {
      if( narg == 3 )
      {
        G4cout << "!!! No input solid type provided. Defaulting to 'box' "
               << G4endl;
      }
      else
      {
        theSolidTypeStr = G4String(vsarg[3]);
      }
        thePosTypeStr = G4String(vsarg[2]);
    }
    thePVTypeStr = G4String(vsarg[1]);
  }
  if( narg > 4 )
  {
    for( G4int ii = 4; ii < narg; ii++ )
    {
      theExtraPars.push_back( vsarg[ii] );
    }
  }

  // Stepping Verbose output class
  G4VSteppingVerbose::SetInstance(new ExDivSteppingVerbose);
  
  // Run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // UserInitialization classes (mandatory)
  ExDivDetectorConstruction* ExDivdetector =
    new ExDivDetectorConstruction(theSolidTypeStr, thePVTypeStr,
                                  thePosTypeStr, theExtraPars );
  runManager->SetUserInitialization(ExDivdetector);
  runManager->SetUserInitialization(new ExDivPhysicsList);
  
  // Set user action classes
  runManager->SetUserInitialization(new ExDivActionInitialization());

#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Initialise G4 kernel
  runManager->Initialize();

  // Get the pointer to the User Interface manager 
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  // Define (G)UI terminal for interactive mode
  if(argc>=1)
  { 
    G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);      
#else
    session = new G4UIterminal();
#endif    

    UI->ApplyCommand("/control/execute vis.mac");    
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

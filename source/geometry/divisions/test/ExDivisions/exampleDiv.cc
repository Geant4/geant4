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
// $Id: exampleDiv.cc,v 1.1 2003-11-19 18:00:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExDivDetectorConstruction.hh"
#include "ExDivPhysicsList.hh"
#include "ExDivPrimaryGeneratorAction.hh"
#include "ExDivRunAction.hh"
#include "ExDivEventAction.hh"
#include "ExDivSteppingAction.hh"
#include "ExDivSteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "ExDivVisManager.hh"
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
    G4cout << "!!! No input solid type provided. Defaulting to 'box' "
           << G4endl;
  }
  else
  {
    if( narg == 2 )
    {
      G4cout << "!!! No input solid type provided. Defaulting to 'box' "
             << G4endl;
    }
    else
    {
      theSolidTypeStr = G4String(vsarg[2]);
    }
    thePVTypeStr = G4String(vsarg[1]);
  }
  if( narg > 3 )
  {
    for( G4int ii = 3; ii < narg; ii++ )
    {
      theExtraPars.push_back( vsarg[ii] );
    }
  }

  // Stepping Verbose output class
  G4VSteppingVerbose::SetInstance(new ExDivSteppingVerbose);
  
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  ExDivDetectorConstruction* ExDivdetector =
    new ExDivDetectorConstruction(theSolidTypeStr, thePVTypeStr, theExtraPars );
  runManager->SetUserInitialization(ExDivdetector);
  runManager->SetUserInitialization(new ExDivPhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new ExDivVisManager;
  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new ExDivPrimaryGeneratorAction(ExDivdetector));
  runManager->SetUserAction(new ExDivRunAction);
  runManager->SetUserAction(new ExDivEventAction);
  runManager->SetUserAction(new ExDivSteppingAction);

  // Initialize G4 kernel
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

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
// This is the *basic* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// To obtain the full version visit the pages: http://sites.google.com/site/hadrontherapy/

// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Main Authors:
//
// G.A.P. Cirrone(a)°, G.Cuttone(a), S.E.Mazzaglia(a), F.Romano(a)
// 
// Contributor authors:
// P.Kaitaniemi(d), A.Heikkinen(d), G.Danielsen (d)
//
// Past authors:
// F.Di Rosa(a), S.Guatelli(c), A.Lechner(e), M.G.Pia(b), G.Russo(a), M.Russo(a), 
//
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
//
// (b) INFN Section of Genova, Italy
// 
// (c) University of Wallongong, Australia
//
// (d) Helsinki Institute of Physics, Helsinki, Finland
//
// (e) CERN, (CH)
//
//  *Corresponding author, email to cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyMatrix.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "globals.hh"
#include "HadrontherapySteppingAction.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyGeometryController.hh"
#include "HadrontherapyGeometryMessenger.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "G4ScoringManager.hh"
#include "IAEAScoreWriter.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc ,char ** argv)
{
  // Set the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
        
  // Only if an initial random seed is needed
  // G4int seed = time(0);
  // CLHEP::HepRandom::setTheSeed(seed);
	
  G4RunManager* runManager = new G4RunManager;
  // Geometry controller is responsible for instantiating the
  // geometries. All geometry specific setup tasks are now in class
  // HadrontherapyGeometryController.
  HadrontherapyGeometryController *geometryController = new HadrontherapyGeometryController();
	
  // Connect the geometry controller to the G4 user interface
  HadrontherapyGeometryMessenger *geometryMessenger = new HadrontherapyGeometryMessenger(geometryController);
	
  G4ScoringManager *scoringManager = G4ScoringManager::GetScoringManager();
  scoringManager->SetVerboseLevel(1);
  scoringManager->SetScoreWriter(new IAEAScoreWriter());
	
  // Initialize the default Hadrontherapy geometry
  geometryController->SetGeometry("default");
	
  // Initialize command based scoring
  G4ScoringManager::GetScoringManager();
	
  // Initialize the physics 
  runManager -> SetUserInitialization(new HadrontherapyPhysicsList());
	
  // Initialize the primary particles
  HadrontherapyPrimaryGeneratorAction *pPrimaryGenerator = new HadrontherapyPrimaryGeneratorAction();
  runManager -> SetUserAction(pPrimaryGenerator);
	
  // Optional UserActions: run, event, stepping
  HadrontherapyRunAction* pRunAction = new HadrontherapyRunAction();
  runManager -> SetUserAction(pRunAction);
	
  HadrontherapyEventAction* pEventAction = new HadrontherapyEventAction();
  runManager -> SetUserAction(pEventAction);
	
  HadrontherapySteppingAction* steppingAction = new HadrontherapySteppingAction(pRunAction); 
  runManager -> SetUserAction(steppingAction);    
	
  // Interaction data: stopping powers
  HadrontherapyInteractionParameters* pInteraction = new HadrontherapyInteractionParameters(true);
	
#ifdef G4VIS_USE
  // Visualization manager
  G4cout << "  SONO DENTRO G4VIS_USE  " << G4endl;
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize();
#endif 
	
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (argc!=1)   // batch mode
    {
      
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);    
     
    }
  else
    {  // interactive mode : define UI session

       
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute defaultMacro.mac");  

#endif
      if (ui->IsGUI())
	UImanager->ApplyCommand("/control/execute macro/GUIPersonalisation.mac");
      ui->SessionStart();
      delete ui;
#endif 
    }

  // Job termination
    // Store dose & fluence data to ASCII & ROOT files 
    if ( HadrontherapyMatrix * pMatrix = HadrontherapyMatrix::GetInstance() )
    {
	pMatrix -> TotalEnergyDeposit(); 
	pMatrix -> StoreDoseFluenceAscii();
#ifdef G4ANALYSIS_USE_ROOT
        pMatrix -> StoreDoseFluenceRoot();
#endif
    }

#ifdef G4ANALYSIS_USE_ROOT
  HadrontherapyAnalysisManager::GetInstance() -> flush();     // Finalize the root file 
#endif


#ifdef G4VIS_USE
  delete visManager;
#endif                


  delete geometryMessenger;
  delete geometryController;
  delete pInteraction; 
  delete runManager;
  return 0;
  
}

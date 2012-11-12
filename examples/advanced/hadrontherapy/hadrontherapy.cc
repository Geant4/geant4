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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
// 
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Main Authors:
//
// R. Calcagno(a), G.A.P. Cirrone(a)*, G.Cuttone(a), F.Romano(a,b), A.Varisano(a)
// 
// Past authors:
// F.Di Rosa(a), S.Guatelli(d), A.Lechner(f), S.E.Mazzaglia(a),  M.G.Pia(c), G.Russo(a), M.Russo(a),
// P.Kaitaniemi(e), A.Heikkinen(e), G.Danielsen (e) 
//
// (a) Laboratori Nazionali del Sud
//     of the INFN, Catania, Italy
//
// (b) Centro Studi e Ricerche e Museo Storico della Fisica E.Fermi, Roma, Italy
// 
// (c) INFN Section of Genova, Italy
// 
// (d) University of Wallongong, Australia
//
// (e) Helsinki Institute of Physics, Helsinki, Finland
//
// (f) CERN, (CH)
//
//  *Corresponding author, email to cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyMatrix.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4UImessenger.hh"
#include "globals.hh"
#include "HadrontherapySteppingAction.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyGeometryController.hh"
#include "HadrontherapyGeometryMessenger.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "G4ScoringManager.hh"


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
        
   G4RunManager* runManager = new G4RunManager;
  // Geometry controller is responsible for instantiating the
  // geometries. All geometry specific setup tasks are now in class
  // HadrontherapyGeometryController.
  HadrontherapyGeometryController *geometryController = new HadrontherapyGeometryController();
	
  // Connect the geometry controller to the G4 user interface
  HadrontherapyGeometryMessenger *geometryMessenger = new HadrontherapyGeometryMessenger(geometryController);
	
  G4ScoringManager *scoringManager = G4ScoringManager::GetScoringManager();
  scoringManager->SetVerboseLevel(1);
  
	
  // Initialize the default Hadrontherapy geometry
  geometryController->SetGeometry("default");
	
  // Initialize command based scoring
  G4ScoringManager::GetScoringManager();
	
  // Initialize the physics 
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = 0;
  G4String physName = ""; 

  // Physics List name defined via environment variable
  char* path = getenv("PHYSLIST");
  if (path) { physName = G4String(path); }

  if(physName != "" && factory.IsReferencePhysList(physName))
    {
      phys = factory.GetReferencePhysList(physName);
    } 

  if(!phys) { phys = new HadrontherapyPhysicsList(); }

  runManager->SetUserInitialization(phys);
  
  //  runManager -> SetUserInitialization(new HadrontherapyPhysicsList());

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
	
  // Initialize analysis
  HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::GetInstance();
#ifdef G4ANALYSIS_USE_ROOT
  analysis -> book();
#endif

#ifdef G4VIS_USE
  // Visualization manager
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
      if(factory.IsReferencePhysList(physName)) 
	{
	  UImanager->ApplyCommand("/control/execute defaultMacroWithReferencePhysicsList.mac");
	}
      else
	{     
	  UImanager->ApplyCommand("/control/execute defaultMacro.mac");  
	}
      
#endif
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
  if (analysis -> IsTheTFile()) analysis -> flush();     // Finalize & write the root file 
#endif


#ifdef G4VIS_USE
  delete visManager;
#endif                


  delete geometryMessenger;
  delete geometryController;
  delete pInteraction; 
  delete runManager;
  delete analysis;
  return 0;
  
}

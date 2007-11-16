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
   
      std::cout << "Usage: load_gdml <filename>" << std::endl;
      return 0;
   }

   G4GDMLParser parser;

   if (!parser.Read(argv[1])) return 0;

   G4RunManager* runManager = new G4RunManager;
   G4VisManager* visManager = new G4VisExecutive;

   runManager->SetUserInitialization(new DetectorConstruction(parser.GetWorldVolume()));
   runManager->SetUserInitialization(new PhysicsList);
   runManager->SetUserAction(new PrimaryGeneratorAction);

   runManager->Initialize();
   visManager->Initialize();

   G4UImanager* UI = G4UImanager::GetUIpointer();

   G4UIsession * session = new G4UIterminal(new G4UItcsh);

   UI->ApplyCommand("/control/execute vis.mac"); 
   UI->ApplyCommand("/run/verbose 0");
   UI->ApplyCommand("/event/verbose 0");
   UI->ApplyCommand("/tracking/verbose 1");
   
   int numberOfEvent = 1;
   runManager->BeamOn(numberOfEvent);

   session->SessionStart();
   
   delete visManager;
   delete runManager;

   return 0;
}

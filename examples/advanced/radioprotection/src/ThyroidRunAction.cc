#include "ThyroidRunAction.hh"
#include "ThyroidEventAction.hh"
#include "ThyroidAnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "ThyroidDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Timer.hh"
ThyroidRunAction::ThyroidRunAction()
{
  //SDname=SDNAME;
 

 
}

ThyroidRunAction::~ThyroidRunAction()
{  
  
}
void ThyroidRunAction::BeginOfRunAction(const G4Run* aRun)
{
 
  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  
 
   ThyroidAnalysisManager* analysis = ThyroidAnalysisManager::getInstance();
   analysis->book();


}




void ThyroidRunAction::EndOfRunAction(const G4Run* aRun)
{
  ThyroidAnalysisManager* analysis = ThyroidAnalysisManager::getInstance();

 



   G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  
   
      analysis->finish();

}





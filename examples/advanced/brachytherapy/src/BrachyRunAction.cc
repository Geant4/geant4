#include "BrachyRunAction.hh"
#include "BrachyEventAction.hh"
#include "BrachyAnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "BrachyDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Timer.hh"
BrachyRunAction::BrachyRunAction(G4String &SDNAME)
{
  SDname=SDNAME;
  //pDetector= (const BrachyDetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  //pEvent=(const BrachyEventAction*)G4RunManager::GetRunManager()->GetUserEventAction();


 
}

BrachyRunAction::~BrachyRunAction()
{  
  
}
void BrachyRunAction::BeginOfRunAction(const G4Run* aRun)
{
 
  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  
 
   BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
   analysis->book();


}




void BrachyRunAction::EndOfRunAction(const G4Run* aRun)
{
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();

 



   G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  
   
      analysis->finish();

}





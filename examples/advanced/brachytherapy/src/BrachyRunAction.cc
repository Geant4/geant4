#include "BrachyRunAction.hh"
#include "BrachyEventAction.hh"
#include "BrachyAnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "BrachyDetectorConstruction.hh"
#include "BrachyRunMessenger.hh"
#include "G4SDManager.hh"
#include "G4Timer.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyFactoryI.hh"
#include "BrachyFactory.hh"
#include "BrachyRunAction.hh"
BrachyRunAction::BrachyRunAction(G4String &SDNAME)
{
  SDname=SDNAME;
  pDet=new BrachyDetectorConstruction(SDname);
  pRun= new BrachyRunMessenger(this);
}

BrachyRunAction::~BrachyRunAction()
{ delete pDet; 
  delete pDetector;
  delete pRun;
}
void BrachyRunAction::BeginOfRunAction(const G4Run* aRun)
{ BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
   analysis->book();
  
   G4RunManager*  pRunManager=G4RunManager::GetRunManager() ;

    if(pRunManager)
     { switch(a)
       { case 1:
              factory=new BrachyFactoryI;
       break;

       default:   
          factory=new BrachyFactoryIr; 
       }
      
       G4VUserPrimaryGeneratorAction* irPrimary=factory->CreatePrimaryGeneratorAction();
      
  
      if(irPrimary)
        { pRunManager->SetUserAction (irPrimary);}
        
     }
}

void BrachyRunAction::SelectEnergy(G4int choice)
{a=choice;
if (a==1)factory=new BrachyFactoryI;
  else factory=new BrachyFactoryIr; 
 
}

void BrachyRunAction::EndOfRunAction(const G4Run* aRun)
{
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();

 



   G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  
      
      analysis->finish();
      delete factory;

      
}





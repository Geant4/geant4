#include "BrachyRunAction.hh"
#include "BrachyEventAction.hh"
#include "BrachyAnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "BrachyDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Timer.hh"
BrachyRunAction::BrachyRunAction(G4String &SDNAME)
{
  SDname=SDNAME;
 pDetector= (const BrachyDetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
 pEvent=(const BrachyEventAction*)G4RunManager::GetRunManager()->GetUserEventAction();

 m_NumVoxelX=pDetector->GetNumVoxelX();
 m_NumVoxelZ=pDetector->GetNumVoxelZ();

 VoxelWidth_X = pDetector->VoxelWidth_X();
 VoxelWidth_Z = pDetector->VoxelWidth_Z();
 
}

BrachyRunAction::~BrachyRunAction()
{  if(pVoxel)delete[]pVoxel;
   delete pDetector;
   delete pEvent;
  
}
void BrachyRunAction::BeginOfRunAction(const G4Run* aRun)
{
 pVoxel=new G4float[ m_NumVoxelX * m_NumVoxelZ];
    for(G4int j=0;j<m_NumVoxelX*m_NumVoxelZ;j++)
                            pVoxel[j] =0.0F;
  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  /* if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
  */
// Book histograms and ntuples
 
   BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
   analysis->book();


}




void BrachyRunAction::EndOfRunAction(const G4Run* aRun)
{
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();

 


 //if (G4VVisManager::GetConcreteInstance()) {
 //  G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
 //        }
   G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
 
  if(pVoxel)
    {

      for(G4int k=0;k<m_NumVoxelZ;k++)
        { 
            z = (-m_NumVoxelZ+1+2*k)*VoxelWidth_Z/2;
                 for(G4int i=0;i<m_NumVoxelX;i++)
		    {  j = i+k*m_NumVoxelX;

		 
                     x = (- m_NumVoxelX +1+2*i)*VoxelWidth_X/2;
		     pVoxel[j]=pEvent->GetEnergy(j); 
		  
		   
		   
       		    if(fabs(x) > 1*mm || fabs(z) > 1*mm)	
		      { 
			 
			 if(pVoxel[j]!=0)

			   {  
                           analysis->analyse(x,z,pVoxel[j]);	
                           analysis->hist(x,z,pVoxel[j]);	
			
                                    		       
                        }
                                           
		 
		      }
		 }
      }
    
    }
    


      analysis->finish();

}





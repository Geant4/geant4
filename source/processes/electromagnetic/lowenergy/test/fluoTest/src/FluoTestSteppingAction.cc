//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestSteppingAction.hh"
#include "FluoTestDetectorConstruction.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "G4ios.hh"
#ifdef G4ANALYSIS_USE
#include "FluoTestAnalysisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#ifdef G4ANALYSIS_USE
FluoTestSteppingAction::FluoTestSteppingAction(
 FluoTestAnalysisManager *aMgr )
  :
    analysisManager(aMgr)
{
  
 }

#else

FluoTestSteppingAction::FluoTestSteppingAction()
{ }
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSteppingAction::~FluoTestSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4double gammaAtTheDetPre=0;
  //G4double gammaTheta = 0;
  //G4double gammatheta = 0;
  // G4double gammaPhi = 0;
  // G4double gammaphi = 0;
  G4double protonsAtTheDetPre=0;
  G4double gammaLeavingSample=0;
  G4double eleLeavingSample=0;
   G4double protonsLeavSam=0;
  
  G4double gammaBornInSample=0;
  
  G4double eleBornInSample=0;
  //G4double otherParticlesAtDetPre=0;
  //counters for data files  
  
  
  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Sample"){
    
    if(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ) 
      { 
	if ((aStep->GetTrack()->GetDynamicParticle()
	     ->GetDefinition()-> GetParticleName()) == "gamma" ) 
	  { //gammaPhi = aStep->GetTrack()->GetMomentumDirection().phi();
	    // gammaTheta = aStep->GetTrack()->GetMomentumDirection().theta();
	 
	  gammaLeavingSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
	  
#ifdef G4ANALYSIS_USE
	  analysisManager->InsGamLeavSam(gammaLeavingSample/keV);  
	  // analysisManager->InsGamLS(gammaTheta);
	  // analysisManager->InsGamLSP(gammaPhi);
#endif 
	  
	  }
      }
    
  }
  
  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Sample"){
    
    if(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ) 
      { 
	if ((aStep->GetTrack()->GetDynamicParticle()
	     ->GetDefinition()-> GetParticleName()) == "e-" ) 
	  {
	    eleLeavingSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    
#ifdef G4ANALYSIS_USE
	    analysisManager->InsEleLeavSam(eleLeavingSample/keV);  
#endif
	  }
	else if ((aStep->GetTrack()->GetDynamicParticle()
	     ->GetDefinition()-> GetParticleName()) == "proton" )
 {
	    protonsLeavSam = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    
#ifdef G4ANALYSIS_USE
	    analysisManager->InsProtLeavSam(protonsLeavSam/keV);  
#endif
	  }

      }
  }
  
  
  if((aStep->GetTrack()->GetDynamicParticle()
      ->GetDefinition()-> GetParticleName()) == "gamma" )
    
    {if(1== (aStep->GetTrack()->GetCurrentStepNumber()))
      
      {if(0 != aStep->GetTrack()->GetParentID())
	
	{if(aStep->GetTrack()->GetVolume()->GetName() == "Sample")
	  {
	    gammaBornInSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    
	     
#ifdef G4ANALYSIS_USE
	    analysisManager->InsGamBornSample(gammaBornInSample/keV);  
#endif 
	    
	  }
	}
      }
    }
  
  if((aStep->GetTrack()->GetDynamicParticle()
       ->GetDefinition()-> GetParticleName()) == "e-" )
    
    {if(1== (aStep->GetTrack()->GetCurrentStepNumber()))
      
       {if(0 != aStep->GetTrack()->GetParentID())
	 
	 {if(aStep->GetTrack()->GetVolume()->GetName() == "Sample")
	   {
	     eleBornInSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    
#ifdef G4ANALYSIS_USE
	     analysisManager->InsEleBornSample(eleBornInSample/keV);  
#endif 
	     
	   }
	 }
       }
    }
  
  if(aStep->GetTrack()->GetNextVolume()){
    
    if(aStep->GetTrack()->GetVolume()->GetName() == "World"){
      
      if(aStep->GetTrack()->GetNextVolume()->GetName() == "HPGeDetector")
	
	{ 
	  if ((aStep->GetTrack()->GetDynamicParticle()
	       ->GetDefinition()-> GetParticleName()) == "gamma" ) 
	    {//gammaphi = aStep->GetTrack()->GetMomentumDirection().phi();
	    //gammatheta = aStep->GetTrack()->GetMomentumDirection().theta();
	  
	    gammaAtTheDetPre = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    
#ifdef G4ANALYSIS_USE
	    analysisManager->InsGamDetPre(gammaAtTheDetPre/keV);  
#endif 

	    }
	  else if ((aStep->GetTrack()->GetDynamicParticle()
		    ->GetDefinition()-> GetParticleName()) == "proton" ) 
	    {
	      protonsAtTheDetPre = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    
#ifdef G4ANALYSIS_USE
	    analysisManager->InsProtDetPre(protonsAtTheDetPre/keV);  
#endif 

	    }
      }
    }
 }
  
  /*
     {if(1== (aStep->GetTrack()->GetCurrentStepNumber()))
       
       {if(0 != aStep->GetTrack()->GetParentID())
	 
	 {if(aStep->GetTrack()->GetVolume()->GetName() == "Sample")
	   {
	     G4cout<<"generato un positrone"<<G4endl;
	     G4cout<<"la sua energia e' "
		   <<(aStep->GetPreStepPoint()->GetKineticEnergy())/keV
		   <<" keV"<<G4endl;
	     G4cout<<"il processo che lo ha generato e' "
		   <<aStep->GetTrack()->GetCreatorProcess()->
	       GetProcessName()<<G4endl;
	       
	   }
	 }
       }
     }
    */



    /*
  if(aStep->GetTrack()->GetNextVolume()){
 
  if(aStep->GetTrack()->GetVolume()->GetName() == "World"){
     
    if(aStep->GetTrack()->GetNextVolume()->GetName() == "HPGeDetector")
      
      { 
	if ((aStep->GetTrack()->GetDynamicParticle()
	     ->GetDefinition()-> GetParticleName()) == "gamma" ) 
	  {
	    gammaAtTheDetPre = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    
#ifdef G4ANALYSIS_USE
	    analysisManager->InsGamDetPre(gammaAtTheDetPre/keV);  
#endif; 
	  }
      }
  }
  }*/
     /*
    if ((aStep->GetTrack()->GetDynamicParticle()
       ->GetDefinition()-> GetParticleName()) == "gamma" ) 
      
    {
      
    if ("Sample"== (aStep->GetTrack()->GetVolume()->GetName()))
      
      { if(1== (aStep->GetTrack()->GetCurrentStepNumber()))
	
	{if(aStep->GetTrack()->GetCreatorProcess()->
	    GetProcessName()!="Transportation")
	  {gammaBornInSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
#ifdef G4ANALYSIS_USE
          analysisManager->InsGamBornSample(gammaBornInSample/keV);  
#endif; 
}
}
}
  //  if ("Sample"!= (aStep->GetTrack()->GetNextVolume()->GetName()))
	//{if(all indietro)
	  //{ gammaLSBackw = (aStep->GetPreStepPoint()->GetKineticEnergy());
	 // #ifdef G4ANALYSIS_USE
         // analysisManager->InsGamLSBackw(gammaLSBackw/keV);  
	 // #endif; 
	//  }
      
      {
	gammaLeavingSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
#ifdef G4ANALYSIS_USE
	analysisManager->InsGamLeavSam(gammaLeavingSample/keV);  
#endif; 
      }
    
  
    if((aStep->GetTrack()->GetVolume()->GetName()!="HPGeDetector")
       &&
       (aStep->GetTrack()->GetNextVolume()->GetName() == "HPGeDetector" ))
      {gammaAtTheDetPre =(aStep->GetPreStepPoint()->GetKineticEnergy());
      
#ifdef G4ANALYSIS_USE
      analysisManager->InsGamDetPre(gammaAtTheDetPre/keV);  
#endif; 	 
      gammaAtTheDetPost =(aStep->GetPostStepPoint()->GetKineticEnergy());
      
#ifdef G4ANALYSIS_USE
      analysisManager->InsGamDetPost(gammaAtTheDetPost/keV);  
#endif; 	
      }
  
    }
  
  if((aStep->GetTrack()->GetDynamicParticle()
      ->GetDefinition()-> GetParticleName()) != "gamma" ) 
    {
      if((aStep->GetTrack()->GetVolume()->GetName()!="HPGeDetector")
	 &&
	 (aStep->GetTrack()->GetNextVolume()->GetName() == "HPGeDetector" ))
	{otherParticlesAtDetPre =(aStep->GetPreStepPoint()->GetKineticEnergy());
#ifdef G4ANALYSIS_USE
	analysisManager->InsOtherPart(otherParticlesAtDetPre/keV);  
#endif; 
	}
    }
  
*/


}	 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

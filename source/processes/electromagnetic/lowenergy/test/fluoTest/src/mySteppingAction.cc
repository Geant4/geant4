//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "mySteppingAction.hh"
#include "myDetectorConstruction.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#ifdef G4ANALYSIS_USE
#include "myAnalysisManager.hh"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#ifdef G4ANALYSIS_USE
mySteppingAction::mySteppingAction(myDetectorConstruction* DET,
 myAnalysisManager *aMgr )
  :detector(DET),analysisManager(aMgr)
{ }

#else

mySteppingAction::mySteppingAction()
{ }
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

mySteppingAction::~mySteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void mySteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4double gammaKenergy=0;
  G4double ioniEnergy=0;
  G4double photoEnergy=0;
  G4double bremEnergy=0;
  G4double comptEnergy=0;
  G4double convEnergy=0;
  G4double raylEnergy=0;
  G4double gammaOutEnergy=0;
  G4double elecKenergy=0;
  G4double eioniEnergy=0;
  G4double ephotoEnergy=0;
  G4double ebremEnergy=0;
  G4double ecomptEnergy=0;
  G4double econvEnergy=0;
  G4double eraylEnergy=0;
  G4double elecOutEnergy=0;

   if(0 != aStep->GetTrack()->GetParentID())
    
     {  if(1== (aStep->GetTrack()->GetCurrentStepNumber()))
	 
       { if((aStep->GetTrack()->GetDynamicParticle()
                   ->GetDefinition()-> GetParticleName()) == "gamma" )
	 { gammaKenergy= (aStep->GetPreStepPoint()->GetKineticEnergy());
#ifdef G4ANALYSIS_USE
          analysisManager->InsertKEnergy(gammaKenergy/keV);  
#endif; 
	 
        if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnergyIoni")
	  {
	  ioniEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
#ifdef G4ANALYSIS_USE
          analysisManager->InsertIoniEnergy(ioniEnergy/keV);  
#endif; 
	  }
     if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnPhotoElec")
	  {photoEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InsertPhotoEnergy(photoEnergy/keV);  
#endif; 
	  } 
      if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnBrem")
	  {bremEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	  G4cout<<"creato un gamma per brem"
		<<G4endl;
#ifdef G4ANALYSIS_USE
          analysisManager->InsertBremEnergy(bremEnergy/keV);  
#endif; 
	  } 
 if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnCompton")
	  {comptEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InsertComptEnergy(comptEnergy/keV);  
#endif; 
	  } 
	  if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnConversion")
	  {convEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InsertConvEnergy(convEnergy/keV);  
#endif; 
	  }
 
      if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnRayleigh")
	  {raylEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InsertRaylEnergy(raylEnergy/keV);  
#endif; 
	  }
	
	if((aStep->GetPreStepPoint()->GetPhysicalVolume()->
	GetName()=="Sample")
	&&
	(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ))
	 	 
	{        
	  gammaOutEnergy= (aStep->GetPreStepPoint()->GetKineticEnergy());
         
#ifdef G4ANALYSIS_USE
          analysisManager->InsertOutEnergy(gammaOutEnergy/keV);  
#endif
	}
      }  
	  
    if((aStep->GetTrack()->GetDynamicParticle()
                   ->GetDefinition()-> GetParticleName()) == "e-" )
	 { elecKenergy= (aStep->GetPreStepPoint()->GetKineticEnergy());
#ifdef G4ANALYSIS_USE
          analysisManager->InserteKEnergy(elecKenergy/keV);  
#endif; 
	 
        if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnergyIoni")
	  {
	  eioniEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
#ifdef G4ANALYSIS_USE
          analysisManager->InserteIoniEnergy(eioniEnergy/keV);  
#endif; 
	  }
     if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnPhotoElec")
	  {ephotoEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InsertePhotoEnergy(ephotoEnergy/keV);  
#endif; 
	  } 
      if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnBrem")
	  {ebremEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InserteBremEnergy(ebremEnergy/keV);  
#endif; 
	  } 
 if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnCompton")
	  {ecomptEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InserteComptEnergy(ecomptEnergy/keV);  
#endif; 
	  } 
	  if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnConversion")
	  {econvEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InserteConvEnergy(econvEnergy/keV);  
#endif; 
	  }
 
      if (aStep->GetTrack()->GetCreatorProcess()->
		   GetProcessName()=="LowEnRayleigh")
	  {eraylEnergy=(aStep->GetPreStepPoint()->GetKineticEnergy());
	
#ifdef G4ANALYSIS_USE
          analysisManager->InserteRaylEnergy(eraylEnergy/keV);  
#endif; 
	  }
	
	if((aStep->GetPreStepPoint()->GetPhysicalVolume()->
	GetName()=="Sample")
	&&
	(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ))
	 	 
	{        
	  elecOutEnergy= (aStep->GetPreStepPoint()->GetKineticEnergy());
         
#ifdef G4ANALYSIS_USE
          analysisManager->InserteOutEnergy(elecOutEnergy/keV);  
#endif
	
        }
	 }   
       }
     }
}

  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

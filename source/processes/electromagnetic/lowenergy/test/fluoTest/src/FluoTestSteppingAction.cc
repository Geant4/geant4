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
				   //FluoTestDetectorConstruction* DET,
 FluoTestAnalysisManager *aMgr )
  :
  //detector(DET),
    analysisManager(aMgr)
{
  nElectrons = 0;
  nDepElec = 0;
  nElecCreated = 0;
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
  //G4double gammaAtTheDetPre=0;
  //G4double gammaAtTheDetPost=0;
  G4double gammaLeavingSample=0;
  //G4double eleLeavingSample=0;
  // G4double gammaLSBackw=0;
  
  //G4double gammaBornInSample=0;
  
  //  G4double eleBornInSample=0;
  //G4double otherParticlesAtDetPre=0;
  //counters for data files  
  
  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Sample"){
    
    if(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ) 
      { 
	if ((aStep->GetTrack()->GetDynamicParticle()
	     ->GetDefinition()-> GetParticleName()) == "gamma" ) 
	  {
	    gammaLeavingSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    
#ifdef G4ANALYSIS_USE
	    analysisManager->InsGamLeavSam(gammaLeavingSample/keV);  
#endif; 
	  }
      }
  }

  if(aStep->GetTrack()->GetNextVolume()){
    
    if(aStep->GetTrack()->GetVolume()->GetName() == "World"){
      
      if(aStep->GetTrack()->GetNextVolume()->GetName() == "HPGeDetector")
	{ 
	  if ((aStep->GetTrack()->GetDynamicParticle()
	       ->GetDefinition()-> GetParticleName()) == "e-" )
	    {
	      
	      G4cout <<"un elettrone e' entrato nel detector "<<G4endl;
	 
	      G4cout<<" ha energia "<<aStep->GetPreStepPoint()->GetKineticEnergy()
		    <<G4endl;
	      nElectrons = nElectrons+1;
	      G4cout<<" ne sono entrati "<< nElectrons<<G4endl;
	    }
	}
    }
  }

  if(aStep->GetTrack()->GetVolume()->GetName() == "HPGeDetector")
    
    { 
      if ((aStep->GetTrack()->GetDynamicParticle()
	   ->GetDefinition()-> GetParticleName()) == "e-" )
	if(1== (aStep->GetTrack()->GetCurrentStepNumber()))
	  {
	    G4cout<<"questo e' nato nel detector " <<G4endl;
	    nElecCreated = nElecCreated+1;
	    nDepElec = nDepElec+1;
	    G4cout<<"ne sono stati creati "<<nElecCreated<<G4endl;
	  }
    }

  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="HPGeDetector")
    {
      
      if(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ) 
	
      {
	if ((aStep->GetTrack()->GetDynamicParticle()
	     ->GetDefinition()-> GetParticleName()) == "e-" )
	  {
	    G4cout<<"un elettrone e' uscito dal detector"<<G4endl;
	    nDepElec=nDepElec-1;
	    G4cout<<"finora sono morti nel detector "<<nDepElec
		  <<" elettroni"<<G4endl;
	  }
      }
    }
}

#ifdef  G4ANALYSIS_USE

#include "XrayFluoAnalysisManager.hh"
#include "XrayFluoAnalysisMessenger.hh"
#include "G4Step.hh"

XrayFluoAnalysisManager* XrayFluoAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::XrayFluoAnalysisManager()
{
  analysisMessenger = new XrayFluoAnalysisMessenger(this);
  histoManager = createIHistoManager(); 
  factory = Lizard::createNTupleFactory();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::~XrayFluoAnalysisManager() 
{

  delete analysisMessenger; 
  analysisMessenger = 0;
  delete  histoManager;
  histoManager = 0;
  delete factory;
  factory = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager* XrayFluoAnalysisManager::getInstance()
{
  if (instance == 0) instance = new XrayFluoAnalysisManager;
  return instance;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XrayFluoAnalysisManager::book()

{
  histoManager->selectStore("XrayFluo.his");

  // Book histograms
  histoManager->create1D("1","Energy Deposit", 100,0.,10.);
 histoManager->create1D("2","Gamma born in the sample", 100,0.,10.);
 histoManager->create1D("3","Electrons  born in the sample", 100,0.,10.);
 histoManager->create1D("4","Gammas leaving the sample", 100,0.,10.);
 histoManager->create1D("5","Electrons leaving the sample ", 100,0.,10.);
 histoManager->create1D("6","Gammas reaching the detector", 100,0.,10.);
 histoManager->create1D("7","Spectrum of the incident particles", 100,0.,10.);
 histoManager->create1D("8","Protons reaching the detector", 100,0.,10.);
 histoManager->create1D("9","Protons leaving the sample", 100,0.,10.);

 // Book ntuples
  ntuple = factory->createC("XrayFluo.his::1");

 //  Add and bind the attributes to the ntuple
  if ( !( ntuple->addAndBind( "energy", eDep) &&
  	  ntuple->addAndBind( "counts"     , counts   ) ) )

    {
      delete ntuple;
      G4Exception("XrayFluoAnalysisManager::book - Could not addAndBind ntuple");
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::finish()
{

  histoManager->store("1");
  histoManager->store("3");
  histoManager->store("2");
  histoManager->store("4");
  histoManager->store("5");
  histoManager->store("6");
  histoManager->store("7");
  histoManager->store("8");
  histoManager->store("9");


  delete ntuple;
  ntuple = 0;
  G4cout << "Deleted ntuple" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analyseStepping(const G4Step* aStep)
{
  G4double gammaAtTheDetPre=0;
  G4double protonsAtTheDetPre=0;
  G4double gammaLeavingSample=0;
  G4double eleLeavingSample=0;
  G4double protonsLeavSam=0;
  G4double gammaBornInSample=0;
  G4double eleBornInSample=0;
  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Sample"){
    
    if(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ) 
      { 
	if ((aStep->GetTrack()->GetDynamicParticle()
	     ->GetDefinition()-> GetParticleName()) == "gamma" )
	  {
	    gammaLeavingSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    IHistogram1D* h1 = histoManager->retrieveHisto1D("4");
	    h1->fill(gammaLeavingSample/keV);
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
	    
	    IHistogram1D* h2 = histoManager->retrieveHisto1D("5");
	    h2->fill(eleLeavingSample/keV);
	  }
	else if ((aStep->GetTrack()->GetDynamicParticle()
		  ->GetDefinition()-> GetParticleName()) == "proton" )
	  {
	    protonsLeavSam = (aStep->GetPreStepPoint()->GetKineticEnergy());
	     IHistogram1D* h3 = histoManager->retrieveHisto1D("9");
	     h3->fill(protonsLeavSam/keV);
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
	    IHistogram1D* h4 = histoManager->retrieveHisto1D("2");
	    h4->fill(gammaBornInSample);
	    
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
	    IHistogram1D* h5 = histoManager->retrieveHisto1D("3");
	    h5->fill(eleBornInSample);
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
	    {
	      gammaAtTheDetPre = (aStep->GetPreStepPoint()->GetKineticEnergy());
	      IHistogram1D* h6 = histoManager->retrieveHisto1D("6");
	      h6->fill( gammaAtTheDetPre);
	    }
	  else if ((aStep->GetTrack()->GetDynamicParticle()
		    ->GetDefinition()-> GetParticleName()) == "proton" ) 
	    {
	      protonsAtTheDetPre = (aStep->GetPreStepPoint()->GetKineticEnergy());
	      IHistogram1D* h7 = histoManager->retrieveHisto1D("8");
	      h7->fill( protonsAtTheDetPre);
	    }
	}
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analyseEnergyDep(G4double energyDep)
{
  IHistogram1D* h8 = histoManager->retrieveHisto1D("1");
  h8->fill(energyDep/keV);
  counts = 1.;
  ntuple->addRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analysePrimaryGenerator(G4double energy)
{
 IHistogram1D* h9 = histoManager->retrieveHisto1D("7");
	      h9->fill(energy/keV);

}
#endif








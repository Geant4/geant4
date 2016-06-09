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
//
// $Id: XrayFluoAnalysisManager.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 11 Jul 2003 A.Mantero, code cleaning / Plotter-XML addiction
//    Sep 2002 A.Mantero, AIDA3.0 Migration
// 06 Dec 2001 A.Pfeiffer updated for singleton
// 30 Nov 2001 Guy Barrand : migrate to AIDA-2.2.
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------
#ifdef G4ANALYSIS_USE

#include "G4VProcess.hh"
#include <fstream>
#include "G4ios.hh"
#include "XrayFluoAnalysisManager.hh"
#include "G4Step.hh"

XrayFluoAnalysisManager* XrayFluoAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::XrayFluoAnalysisManager()
  :outputFileName("xrayfluo"), visPlotter(false), phaseSpaceFlag(false), physicFlag (false), persistencyType("xml"), 
   deletePersistencyFile(true), gunParticleEnergies(0), gunParticleTypes(0), analysisFactory(0), tree(0),histogramFactory(0), plotter(0)
{
  //creating the messenger
  analisysMessenger = new XrayFluoAnalysisMessenger(this);

  // Hooking an AIDA compliant analysis system.
  // creating Analysis factroy, necessary to create/manage analysis
  analysisFactory = AIDA_createAnalysisFactory();

  CreatePersistency(outputFileName,persistencyType);

  if(analysisFactory) {

    
    // Creating the plotter factory

    if (visPlotter){
      plotterFactory = analysisFactory->createPlotterFactory();
    }
  }
  
  G4cout << "XrayFluoAnalysisManager created" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::~XrayFluoAnalysisManager() 
{
  delete histogramFactory;
  histogramFactory=0;

  delete analysisFactory;
  analysisFactory = 0;

  delete plotterFactory;
  plotterFactory=0;

  delete tree;

  if ( gunParticleEnergies ) delete gunParticleEnergies;
  gunParticleEnergies = 0;
  if ( gunParticleTypes ) delete gunParticleTypes;
  gunParticleTypes = 0;

  delete instance;

  G4cout << "XrayFluoAnalysisManager delete" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager* XrayFluoAnalysisManager::getInstance()

{
  if (instance == 0) {instance = new XrayFluoAnalysisManager;}
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::CreatePersistency(G4String fileName,G4String persistencyType,
						G4bool readOnly, G4bool createNew)
{

  if (tree) {
    delete tree;
    tree = 0;
  }
  if (treeDet) {
    delete treeDet;
    treeDet =0;
  }
  if (histogramFactory) {
    delete histogramFactory;
    histogramFactory = 0;
  }
  if (tupleFactory) {
    delete tupleFactory;
    tupleFactory = 0; 
  }

  if (tupleDetFactory) {
    delete tupleDetFactory;
    tupleDetFactory = 0;
  }
  G4String fileNameDet = fileName+"Detector";
  
  AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
  if(treeFactory) {
    if (persistencyType == "hbook") {
      fileName = fileName + ".hbk";
      fileNameDet = fileNameDet + ".hbk";
    }
    else if (persistencyType == "xml"){
      fileName = fileName + ".xml";
      fileNameDet = fileNameDet + ".xml";
    }
    
    if (phaseSpaceFlag) {
      
      tree = treeFactory->create(fileName,persistencyType,readOnly,createNew); // output file
      treeDet = treeFactory->create(fileNameDet,persistencyType,readOnly,createNew); // output file
      if(analysisFactory) {
       	
	tupleDetFactory = analysisFactory->createTupleFactory(*treeDet);  
	tupleFactory = analysisFactory->createTupleFactory(*tree);
	
      }
    } 
    
    // trees to be stored in case of phase-space production
    else {
      
      tree = treeFactory->create(fileName,persistencyType,readOnly,createNew); // output file
      
      if(analysisFactory) {
       	
	histogramFactory = analysisFactory->createHistogramFactory(*tree);  
	
      }
    }
    
    delete treeFactory; // Will not delete the ITree.
  }
}


void XrayFluoAnalysisManager::book()
{
  
  if (phaseSpaceFlag) {
    /*    
    // Book clouds
    
    cloud_1 = histogramFactory->createCloud1D("Gamma Exting Sample","ciao!",-1); 
    cloud_2 = histogramFactory->createCloud1D("Gamma Incident on the detector","ciao!",-1);
    cloud_3 = histogramFactory->createCloud1D("Electrons Exiting the Sample","ciao!",-1);
    beamCloud = histogramFactory->createCloud1D("Incident Radiation Spectrum","ciao!",-1);
    */
    // Book output Tuple

    // Book tuple
    std::vector<std::string> columnNames;
    columnNames.push_back("Particle");
    columnNames.push_back("Energies");
    columnNames.push_back("momentumTheta");
    columnNames.push_back("momentumPhi");
    columnNames.push_back("Processes");

    std::vector<std::string> columnTypes;
    columnTypes.push_back("int");
    columnTypes.push_back("double");
    columnTypes.push_back("double");
    columnTypes.push_back("double");
    columnTypes.push_back("int"); // useful for hbk


    tupleFluo = tupleFactory->create("10", "Total Tuple", columnNames, columnTypes, "");
    assert(tupleFluo);

  }
  
  else {
    // Book histograms
    
    histo_1 = histogramFactory->createHistogram1D("1","Energy Deposit", 500,0.,10.); //20eV def.
    histo_2 = histogramFactory->createHistogram1D("2","Gamma born in the sample", 100,0.,10.);
    histo_3 = histogramFactory->createHistogram1D("3","Electrons  born in the sample", 100,0.,10.);
    histo_4 = histogramFactory->createHistogram1D("4","Gammas leaving the sample", 300,0.,10.);
    histo_5 = histogramFactory->createHistogram1D("5","Electrons leaving the sample ",200000 ,0.,10.0); // .05 eV def.
    histo_6 = histogramFactory->createHistogram1D("6","Gammas reaching the detector", 100,0.,10.);
    histo_7 = histogramFactory->createHistogram1D("7","Spectrum of the incident particles", 100,0.,10.);
    histo_8 = histogramFactory->createHistogram1D("8","Protons reaching the detector", 100,0.,10.);
    histo_9 = histogramFactory->createHistogram1D("9","Protons leaving the sample", 100,0.,10.);
    
    // Debugging-purpose Histos
    
    //histo_10 = histogramFactory->createHistogram1D("10","Photon Origin", 4,0.,3.);  
    //histo_11 = histogramFactory->createHistogram1D("11","Spectrum from LowEnPhotoELectric", 300,0.,10.);
    //histo_12 = histogramFactory->createHistogram1D("12","Spectrum From the other processes (unknown)", 300,0.,10.);
    
    //  IHistogram2D* histo_20 = histogramFactory->
    //create2D("20","Phi, Theta",80 ,-3.14,3.14,80,0.,3.14);  
  }
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::LoadGunData(G4String fileName, G4bool raileighFlag) {
  
  //  gunParticleEnergies = new std::vector<G4double>;
  //  gunParticleTypes = new std::vector<G4String>;

  G4String ext = fileName.substr(fileName.size()-3,fileName.size()-1);
  G4String persistencyType;
  
  if (ext == "xml") {
    
    persistencyType = "xml";
  }
  else if (ext == "hbk") {
    
    persistencyType = "hbook";
  }

  gunParticleEnergies = new std::vector<G4double>;
  gunParticleTypes = new std::vector<G4String>;

  AIDA::ITreeFactory* treeDataFactory = analysisFactory->createTreeFactory();
  AIDA::ITree* treeData = treeDataFactory->create(fileName,persistencyType,true,false); // input file
  AIDA::IManagedObject* mo = treeData->find("10");
  AIDA::ITuple* tupleData = dynamic_cast<AIDA::ITuple*>(mo);
  tupleData->start();

  while (tupleData->next()) {
    if (raileighFlag ^ (!raileighFlag && (tupleData->getInt(4)) ) ) {
      gunParticleEnergies->push_back(tupleData->getDouble(1));
      if (tupleData->getInt(0) == 1 ) gunParticleTypes->push_back("gamma");
      if (tupleData->getInt(0) == 0 ) gunParticleTypes->push_back("e-");
    }

  }
  G4cout << "Maximum mumber of events: "<< gunParticleEnergies->size() <<G4endl;


}

std::vector<G4double>* XrayFluoAnalysisManager::GetEmittedParticleEnergies() {

  return gunParticleEnergies;

}

std::vector<G4String>* XrayFluoAnalysisManager::GetEmittedParticleTypes() {

  return gunParticleTypes;

}


void XrayFluoAnalysisManager::finish()
{

  if (tupleFluo) {ExtractData();};

  if(plotter)  {
    // Wait for the keyboard return to avoid destroying the plotter window too quickly.
    G4cout << "Press <ENTER> to exit" << G4endl;
    G4cin.get();
    plotter->hide(); //hide plotter windows, but doesn't delete plotter
  }
  
  if(tree) {
    tree->commit(); // Write histos in file. 
    tree->close();
  }
  if (treeDet) {
    treeDet->commit();
    treeDet->close();
  }
  
  deletePersistencyFile = false;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::InitializePlotter()
{

  // If speciefied (visPlotter=true),
  // a window for the visulizaton of partial results is created
  if(plotterFactory && visPlotter && deletePersistencyFile) 
    {
      plotter = plotterFactory->create();
      // Set the page title :
      plotter->setParameter("pageTitle","XrayFluo");
    }

  if(plotter && visPlotter) {
    plotter->show(); // shows plotter window
  }
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::PlotCurrentResults()
{
  if(plotter) {
    if (phaseSpaceFlag){
      // Plotting the spectrum of Gamma exiting the sample - cloud1
      AIDA::ICloud1D& cloud = *cloud_1;
      AIDA::IFilter* filterGamma = tupleFactory->createFilter(
                                              " Particle == std::string(\"gamma\") ");
      AIDA::IEvaluator* evaluatorEnergy = tupleFactory->createEvaluator("Energies");
      filterGamma->initialize(*tupleFluo); 
      evaluatorEnergy->initialize(*tupleFluo);
      tupleFluo->project(cloud,*evaluatorEnergy,*filterGamma);  
      
      plotter->currentRegion().plot( cloud, "Exiting Gammas " );
      plotter->refresh();
    }
    
    else{
      // Plotting the Detector Energy Deposit - histo_1
      AIDA::IHistogram1D& histo1p = *histo_1;
      plotter->currentRegion().plot( histo1p, "Detector Energy Deposition" );
      plotter->refresh();
    }
    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool XrayFluoAnalysisManager::GetDeletePersistencyFileFlag()
{
  return deletePersistencyFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::SetPhysicFlag(G4bool val)
{
  physicFlag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





void XrayFluoAnalysisManager::analyseStepping(const G4Step* aStep)
{
  
  if (phaseSpaceFlag){
    
    
    G4String particleType="";
    G4String parentProcess="";
    G4ThreeVector momentum=0;
    G4double particleEnergy=0;
    
    // Select volume from wich the step starts
    if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Sample"){
      // Select volume in wich the step ends
      if (physicFlag ^ (!physicFlag && (aStep->GetTrack()->GetNextVolume()->GetName() == "World" ))) { 
	//
	//
	//	G4cout << "physicFlag: "<< physicFlag << G4endl
	//	       << "NextVolume: "<< aStep->GetTrack()->GetNextVolume()->GetName() << G4endl;
	
	// extracting information needed
	particleType = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
	momentum = aStep->GetTrack()->GetDynamicParticle()->GetMomentum();
	particleEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
	G4int parent;
	if(aStep->GetTrack()->GetCreatorProcess()){
	  parentProcess = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
	  parent = 1;
	}	
	else {
	  parentProcess = "Not Known";
	  parent = 0;
	}
	// filling tuple 
	
	//	G4cout<< particleType << G4endl;
	G4int part = 2 ;
	if (particleType == "gamma") part =1; 
	if (particleType == "e-") part = 0;
	
	tupleFluo->fill(0,part);
	tupleFluo->fill(1,particleEnergy);
	tupleFluo->fill(2,momentum.theta());
	tupleFluo->fill(3,momentum.phi());
	tupleFluo->fill(4,parent); //hacked to be useful for hbk
	
	tupleFluo->addRow();
      }
    }
  }
  
  
  // Normal behaviour, without creation of phase space
  else {
    
  G4double gammaAtTheDetPre=0;
  G4double protonsAtTheDetPre=0;
  G4double gammaLeavingSample=0;
  //G4double gammaLeavingSamplePhi=0;
  //G4double gammaLeavingSampleTheta=0;

  G4double eleLeavingSample=0;
  G4double protonsLeavSam=0;
  G4double gammaBornInSample=0;
  G4double eleBornInSample=0;

  // Filling the histograms with data, passing thru stepping.


  // Select volume from wich the step starts
  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Sample"){
    // Select volume from wich the step ends
    if(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ) { 
      // Select the particle type
      if ((aStep->GetTrack()->GetDynamicParticle()
	   ->GetDefinition()-> GetParticleName()) == "gamma" )
	
	{
	    
	    
	  // Control histos, used for debugging purpose.

	  //  	    if(aStep->GetTrack()->GetCreatorProcess()){
	  //  	      G4String process = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
	  //  	      G4cout << "Il processo di origine e' " << process << G4endl;
	  //  	      if(histo_10) {
	  //  		histo_10->fill(1.);
	  //  		gammaLeavingSample = 
	  //  		  (aStep->GetPreStepPoint()->GetKineticEnergy());
	  //  		if(histo_11) {
	  //  		  histo_11->fill(gammaLeavingSample/keV);
		  
	  //  		}
	  //  	      }    
	  //  	    }
	  //  	    else {
	  //  	      //G4cout << "Sembra che non ci sia un processo d'origine" << G4endl;
	  //  	      if(histo_10) {
	  //  		histo_10->fill(2.);

	  //  		gammaLeavingSample = 
	  //  		  (aStep->GetPreStepPoint()->GetKineticEnergy());
	  //  		if(histo_12) {
	  //  		  histo_12->fill(gammaLeavingSample/keV);

	  //  		}
		
	  //  	      }
	  //  	    }
	    

	  //                  //
	  //  Filling Histos  //
	  //                  // 


	  gammaLeavingSample = 
	    (aStep->GetPreStepPoint()->GetKineticEnergy());
	  if(histo_4) {
	    histo_4->fill(gammaLeavingSample/keV);

	  }

	  // Other debugging purpose histos
	  /*    gammaLeavingSamplePhi = 
		(aStep->GetPreStepPoint()->GetMomentumDirection().phi());
		G4cout << "questo e' Phi: " << gammaLeavingSamplePhi << G4endl;
		gammaLeavingSampleTheta = 
		(aStep->GetPreStepPoint()->GetMomentumDirection().theta());
		G4cout << "questo e' Theta: " << gammaLeavingSampleTheta << G4endl;
		if(histo_20) {
		// G4cout << "histo_20 esiste" << G4endl;
		histo_20->fill(gammaLeavingSamplePhi,gammaLeavingSampleTheta,1.);
		}  */
  

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
	    if(histo_5) {
	      histo_5->fill(eleLeavingSample/keV);
	    }
	  }
	else if ((aStep->GetTrack()->GetDynamicParticle()
		  ->GetDefinition()-> GetParticleName()) == "proton" )
	  {
	    protonsLeavSam = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    if(histo_9) {
	      histo_9->fill(protonsLeavSam/keV);
	    }
	  }
	
      }
  }

  // electrons  from the detector -- Debugging
  
//   if((aStep->GetTrack()->GetDynamicParticle()
//       ->GetDefinition()-> GetParticleName()) == "e-" ){
    
//     if(1== (aStep->GetTrack()->GetCurrentStepNumber())){ 
      
//       if(0 != aStep->GetTrack()->GetParentID()){
// 	if(aStep->GetTrack()->GetVolume()->GetName() == "HPGeDetector")
	  
// 	  if(aStep->GetTrack()->GetNextVolume()->GetName() == "World" ){
	    
// 	    if(aStep->GetTrack()->GetCreatorProcess()){
//  	      G4String process = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
//  	      G4cout << "Origin Porcess Name:  " << process << G4endl;
//  	      if(histo_10) {
//  		histo_10->fill(1.);
//  		gammaLeavingSample = 
//  		  (aStep->GetPreStepPoint()->GetKineticEnergy());
//  		if(histo_11) {
//  		  histo_11->fill(gammaLeavingSample/keV);
		  
//  		}
//  	      }    
//  	    }
//  	    else {
//  	      G4cout << "No Origin Process Name" << G4endl;
//  	      if(histo_10) {
//  		histo_10->fill(2.);
		
//  		gammaLeavingSample = 
//  		  (aStep->GetPreStepPoint()->GetKineticEnergy());
//  		if(histo_12) {
//  		  histo_12->fill(gammaLeavingSample/keV);
		  
//  		}
		
//  	      }
//  	    }
// 	  }
//       }
      
//     }
//   }



 
  if((aStep->GetTrack()->GetDynamicParticle()
      ->GetDefinition()-> GetParticleName()) == "gamma" )
    
    {if(1== (aStep->GetTrack()->GetCurrentStepNumber()))
      
      {if(0 != aStep->GetTrack()->GetParentID())
	
	{if(aStep->GetTrack()->GetVolume()->GetName() == "Sample")
	  {
	    gammaBornInSample = (aStep->GetPreStepPoint()->GetKineticEnergy());
	    if(histo_2) {
	      histo_2->fill(gammaBornInSample/keV);
	    }
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
	    if(histo_3) {
	      histo_3->fill(eleBornInSample/keV);
	    }
	  }
	}
      }
    }
  
  if(aStep->GetTrack()->GetNextVolume()){
    
    //if(aStep->GetTrack()->GetVolume()->GetName() == "World"){
      
    if(aStep->GetTrack()->GetNextVolume()->GetName() == "HPGeDetector")
	
      { 
	if ((aStep->GetTrack()->GetDynamicParticle()
	     ->GetDefinition()-> GetParticleName()) == "gamma" ) 
	  {

	    gammaAtTheDetPre = 
	      (aStep->GetPreStepPoint()->GetKineticEnergy());
	    if(histo_6) {
	      histo_6->fill( gammaAtTheDetPre/keV);
	    }
	  }
	else if ((aStep->GetTrack()->GetDynamicParticle()
		  ->GetDefinition()-> GetParticleName()) == "proton" ) 
	  {
	    protonsAtTheDetPre = 
	      (aStep->GetPreStepPoint()->GetKineticEnergy());
	    if(histo_8) {
	      histo_8->fill( protonsAtTheDetPre/keV);
	    }
	  }
      }
    //}  close of if(aStep->GetTrack()->GetVolume()->GetName() == "World"){ 
  }
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::ExtractData(){
  
  if (tupleFluo->rows()) {
    
    
    //    AIDA::IFilter* filterGamma = tupleFactory->createFilter(" Particle == std::string(\"gamma\")");
    //    AIDA::IFilter* filterEminus = tupleFactory->createFilter(" Particle == std::string(\"e-\")");
    
    
    AIDA::IFilter* filterAngle = tupleFactory->createFilter("(momentumPhi   >= (220. * (3.1415926/180.) )) && 
                                                             (momentumPhi   <= (230. * (3.1415926/180.) )) && 
                                                             (momentumTheta >= (130. * (3.1415926/180.) )) && 
                                                             (momentumTheta <= (140. * (3.1415926/180.) )) " );
    
    
    //    filterGamma  ->initialize(*tupleFluo); 
    //    filterEminus ->initialize(*tupleFluo);
    filterAngle->initialize(*tupleFluo);
    
    // Create IEvaluator and initialize it to this ITuple
    //    AIDA::IEvaluator* evaluatorEnergy = tupleFactory->createEvaluator("Energies");
    //    evaluatorEnergy->initialize(*tupleFluo);
    
    //    tupleFluo->project(*cloud_1,*evaluatorEnergy,*filterGamma);  
    //    tupleFluo->project(*cloud_2,*evaluatorEnergy,*filterAngle);  
    //    tupleFluo->project(*cloud_3,*evaluatorEnergy,*filterEminus);  
    
    tupleDetFluo = tupleDetFactory->createFiltered("1", *tupleFluo, *filterAngle);
    assert(tupleDetFluo);    
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XrayFluoAnalysisManager::analyseEnergyDep(G4double energyDep)
{

  // Filling of Energy Deposition, called by XrayFluoEventAction
  
  histo_1->fill(energyDep/keV);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analysePrimaryGenerator(G4double energy)
{

  // Filling of energy spectrum histogram of the primary generator

  if (phaseSpaceFlag){
    //    beamCloud->fill(energy/keV);
  }
  else {
    histo_7->fill(energy/keV);
  }
}


// not used -- Created for future development

void XrayFluoAnalysisManager::SetOutputFileName(G4String newName)
{

  outputFileName = newName;
}

void XrayFluoAnalysisManager::SetOutputFileType(G4String newType)
{

  persistencyType = newType;
}

#endif









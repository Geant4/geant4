//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
#include <stdlib.h>
#include "G4VProcess.hh"
#include <fstream>
#include <strstream>
#include "G4ios.hh"
#include "XrayFluoAnalysisManager.hh"
#include "G4Step.hh"

XrayFluoAnalysisManager* XrayFluoAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::XrayFluoAnalysisManager()
  :outputFileName("xrayfluo"), /*visPlotter(false),*/ persistencyType("xml"), 
  deletePersistencyFile(true), analysisFactory(0), tree(0),histogramFactory(0)//, plotter(0)
{
  //creating the messenger
  analisysMessenger = new XrayFluoAnalysisMessenger(this);

  // Hooking an AIDA compliant analysis system.
  // creating Analysis factroy, necessary to create/manage analysis
  analysisFactory = AIDA_createAnalysisFactory();

    // creating  persistency

  CreatePersistency(outputFileName,persistencyType);

//   if(analysisFactory) {


    
//     // Creating the plotter factory
//     plotterFactory = analysisFactory->createPlotterFactory();
//   }


  G4cout << "XrayFluoAnalysisManager created" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::~XrayFluoAnalysisManager() 
{
  delete histogramFactory;
  histogramFactory=0;

  delete analysisFactory;
  analysisFactory = 0;

  //  delete plotterFactory;
  //  plotterFactory=0;

  delete tree;

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

  if (tree) delete tree;
  if (histogramFactory) delete histogramFactory;

    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if(treeFactory) {
      if (persistencyType == "hbook") {
	fileName = fileName + ".hbk";
      }
      else if (persistencyType == "xml"){
	fileName = fileName + ".xml";
      }
      tree = treeFactory->create(fileName,persistencyType,readOnly,createNew); // output file
      
      delete treeFactory; // Will not delete the ITree.
    }
  if(analysisFactory) {


    histogramFactory = analysisFactory->createHistogramFactory(*tree);  

  }
}


void XrayFluoAnalysisManager::book()
{
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
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::finish()
{

//   if(plotter)  {
//     // Wait for the keyboard return to avoid destroying the plotter window too quickly.
//     G4cout << "Press <ENTER> to exit" << G4endl;
//     G4cin.get();
//     plotter->hide(); //hide plotter windows, but doesn't delete plotter
//   }

  if(tree) {
    tree->commit(); // Write histos in file. 
    tree->close();
  }

  deletePersistencyFile = false;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void XrayFluoAnalysisManager::InitializePlotter()
// {

//   // If speciefied (visPlotter=true),
//   // a window for the visulizaton of partial results is created
//   if(plotterFactory && visPlotter && deletePersistencyFile) 
//     {
//       plotter = plotterFactory->create();
//       // Set the page title :
//       plotter->setParameter("pageTitle","XrayFluo");
//     }

//   if(plotter && visPlotter) {
//     plotter->show(); // shows plotter window
//   }
  
// }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void XrayFluoAnalysisManager::PlotCurrentResults()
// {
//   if(plotter) {
//       // Plotting the Detector Energy Deposit - histo_1
//       AIDA::IHistogram1D& histo1p = *histo_1;
//       plotter->currentRegion().plot( histo1p, "Detector Energy Deposition" );
//       plotter->refresh();
//     }
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool XrayFluoAnalysisManager::GetDeletePersistencyFileFlag()
{
  return deletePersistencyFile;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analyseStepping(const G4Step* aStep)
{

  
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

    histo_7->fill(energy/keV);
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











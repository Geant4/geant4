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
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 06 Dec 2001 A.Pfeiffer updated for singleton
// 30 Nov 2001 Guy Barrand : migrate to AIDA-2.2.
// 28 Nov 2001 Elena Guardincerri     Created
// 29 Nov 2002 Muigration to AIDA 3.0 (Alfonso.mantero@ge.infn.it)

// -------------------------------------------------------------------
#include <stdlib.h>
#include "G4VProcess.hh"
#include "g4std/fstream"
#include "G4ios.hh"
#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

#include "XrayFluoAnalysisManager.hh"


#include "G4Step.hh"

XrayFluoAnalysisManager* XrayFluoAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::XrayFluoAnalysisManager()
  :outputFileName("xrayfluo.hbk"),analysisFactory(0), tree(0),histogramFactory(0)
{
  //  pEvent=new XrayFluoEventAction();

  analisysMessenger = new XrayFluoAnalysisMessenger(this);

  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory) {

    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if(treeFactory) {
      // tree = treeFactory->create(); // Tree in memory.
      //porebbe essere qua il memory leak//

      tree = treeFactory->create(outputFileName,"hbook",false,true);

      delete treeFactory; // Will not delete the ITree.
      histogramFactory = analysisFactory->createHistogramFactory(*tree);  
      // tupleFactory = analysisFactory->createTupleFactory(*tree);
    }

    /*
    // Create and a the plotter :
    IPlotterFactory* plotterFactory = 
      analysisFactory->createPlotterFactory(argc,argv);
    if(plotterFactory) {
      plotter = plotterFactory->create();
      if(plotter) {
        // Map the plotter on screen :
	plotter->show();
        // Set the page title :
	plotter->setParameter("pageTitle","XrayFluo");
	// Have two plotting regions (one column, two rows).
      }
      delete plotterFactory;
    }
    */
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

  // delete tupleFactory;
  // tupleFactory=0;

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

void XrayFluoAnalysisManager::book()
{
    // Book histograms

  histo_1 = histogramFactory->createHistogram1D("1","Energy Deposit", 100000,0.,10.); //1eV def.
  histo_2 = histogramFactory->createHistogram1D("2","Gamma born in the sample", 100,0.,10.);
  histo_3 = histogramFactory->createHistogram1D("3","Electrons  born in the sample", 100,0.,10.);
  histo_4 = histogramFactory->createHistogram1D("4","Gammas leaving the sample", 300,0.,10.);
  histo_5 = histogramFactory->createHistogram1D("5","Electrons leaving the sample ",200000 ,0.,10.0); // .05 eV def.
  histo_6 = histogramFactory->createHistogram1D("6","Gammas reaching the detector", 100,0.,10.);
  histo_7 = histogramFactory->createHistogram1D("7","Spectrum of the incident particles", 100,0.,10.);
  histo_8 = histogramFactory->createHistogram1D("8","Protons reaching the detector", 100,0.,10.);
  histo_9 = histogramFactory->createHistogram1D("9","Protons leaving the sample", 100,0.,10.);
  histo_10 = histogramFactory->createHistogram1D("10","Photon Origin", 4,0.,3.);  
  histo_11 = histogramFactory->createHistogram1D("11","Spectrum from LowEnPhotoELectric", 300,0.,10.);
  histo_12 = histogramFactory->createHistogram1D("12","Spectrum From the other processes (unknown)", 300,0.,10.);

  //  IHistogram2D* histo_20 = histogramFactory->
  //create2D("20","Phi, Theta",80 ,-3.14,3.14,80,0.,3.14);  
  
  // Create a tuple :
  // tuple = tupleFactory->create("XrayFluo","XrayFluo","energy counts");
} 
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::finish()
{

  if(tree) {
    tree->commit(); // Write histos and tuple in file. 
    tree->close();
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analyseStepping(const G4Step* aStep)
{

 
 /*
  IHistogram1D*   histo_2 =dynamic_cast<IHistogram1D *> ( tree->find("2") );
  IHistogram1D*   histo_3 =dynamic_cast<IHistogram1D *> ( tree->find("3") );
  IHistogram1D*   histo_4 =dynamic_cast<IHistogram1D *> ( tree->find("4") );
  IHistogram1D*   histo_5 =dynamic_cast<IHistogram1D *> ( tree->find("5") );
  IHistogram1D*   histo_6 =dynamic_cast<IHistogram1D *> ( tree->find("6") );
  IHistogram1D*   histo_8 =dynamic_cast<IHistogram1D *> ( tree->find("8") );
  IHistogram1D*   histo_10 =dynamic_cast<IHistogram1D *> ( tree->find("10") );  
  IHistogram1D*   histo_11 =dynamic_cast<IHistogram1D *> ( tree->find("11") );
  IHistogram1D*   histo_12 =dynamic_cast<IHistogram1D *> ( tree->find("12") );
  //IHistogram2D*   histo_20=dynamic_cast<IHistogram2D *> ( tree->find("20") );
  */
  
  G4double gammaAtTheDetPre=0;
  G4double protonsAtTheDetPre=0;
  G4double gammaLeavingSample=0;
  //G4double gammaLeavingSamplePhi=0;
  //G4double gammaLeavingSampleTheta=0;

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
	    
	    
	    if(aStep->GetTrack()->GetCreatorProcess()){
	      G4String process = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
	      // G4cout << "Il processo di origine e' " << process << G4endl;
	      if(histo_10) {
		histo_10->fill(1.);
		gammaLeavingSample = 
		  (aStep->GetPreStepPoint()->GetKineticEnergy());
		if(histo_11) {
		  histo_11->fill(gammaLeavingSample/keV);
		  
		}
	      }    
	    }
	    else {
	      //G4cout << "Sembra che non ci sia un processo d'origine" << G4endl;
	      if(histo_10) {
		histo_10->fill(2.);

		gammaLeavingSample = 
		  (aStep->GetPreStepPoint()->GetKineticEnergy());
		if(histo_12) {
		  histo_12->fill(gammaLeavingSample/keV);

		}
		
	      }
	    }
	    


	    gammaLeavingSample = 
	      (aStep->GetPreStepPoint()->GetKineticEnergy());
	    if(histo_4) {
	      histo_4->fill(gammaLeavingSample/keV);

	    }

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
	    if(histo_3) {
	      histo_3->fill(protonsLeavSam/keV);
	    }
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
    
    if(aStep->GetTrack()->GetVolume()->GetName() == "World"){
      
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
    }
    }

  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analyseEnergyDep(G4double energyDep)
{
  //histo_1 =dynamic_cast<IHistogram1D *> ( tree->find("1") );
  
  histo_1->fill(energyDep/keV);
  
  /*
    if(tuple) {
    tuple->fill(0,energyDep);
    tuple->fill(1,1);
    tuple->addRow();
    }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analysePrimaryGenerator(G4double energy)
{

  // qua si fa uso spropositato di memoria -- attenzione!!//
  //histo_7 =dynamic_cast<IHistogram1D *> ( tree->find("7") );
    histo_7->fill(energy/keV);
}

void XrayFluoAnalysisManager::SetOutputFileName(G4String newName)
{

  outputFileName = newName;

}













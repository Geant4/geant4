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
// 06 Dec 2001 A.Pfeiffer updated for singleton
// 30 Nov 2001 Guy Barrand : migrate to AIDA-2.2.
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------
#include <stdlib.h>
#include "g4std/fstream"
#include "G4ios.hh"
#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

#include "XrayFluoAnalysisManager.hh"


#include "G4Step.hh"

XrayFluoAnalysisManager* XrayFluoAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::XrayFluoAnalysisManager()
  :analysisFactory(0), tree(0),histogramFactory(0) ,tupleFactory(0)
{
  //  pEvent=new XrayFluoEventAction();


  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory) {

    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if(treeFactory) {
      // tree = treeFactory->create(); // Tree in memory.
      tree = treeFactory->create("XrayFluo.hbk",false,false,"hbook");

      delete treeFactory; // Will not delete the ITree.
     histogramFactory = analysisFactory->createHistogramFactory(*tree);  
     tupleFactory = analysisFactory->createTupleFactory(*tree);
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisManager::~XrayFluoAnalysisManager() 
{

  


 delete histogramFactory;
  histogramFactory=0;

  delete analysisFactory;
  analysisFactory = 0;

  delete tupleFactory;
  tupleFactory=0;
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
  

 
  // Book histograms
  IHistogram1D*   histo_1 = histogramFactory->create1D
    ("1","Energy Deposit", 100,0.,10.);
  IHistogram1D*  histo_2 = histogramFactory->create1D
    ("2","Gamma born in the sample", 100,0.,10.);
  IHistogram1D*   histo_3 = histogramFactory->create1D
    ("3","Electrons  born in the sample", 100,0.,10.);
  IHistogram1D*   histo_4 = histogramFactory->create1D
    ("4","Gammas leaving the sample", 100,0.,10.);
  IHistogram1D*  histo_5 = histogramFactory->create1D
    ("5","Electrons leaving the sample ", 100,0.,10.);
  IHistogram1D*   histo_6 = histogramFactory->create1D
    ("6","Gammas reaching the detector", 100,0.,10.);
  IHistogram1D*   histo_7 = histogramFactory->create1D
    ("7","Spectrum of the incident particles", 100,0.,10.);
  IHistogram1D*   histo_8 = histogramFactory->create1D
    ("8","Protons reaching the detector", 100,0.,10.);
  IHistogram1D*  histo_9 = histogramFactory->create1D
    ("9","Protons leaving the sample", 100,0.,10.);
  IHistogram1D*   histo_10 = histogramFactory->create1D
    ("10","Gammas leaving the sample smeared", 100,0.,10.);  
  
  


 // Create a tuple :
// tuple = tupleFactory->create("XrayFluo","XrayFluo","energy counts");
	
} 




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::finish()
{

  if(tree) {tree->commit(); // Write histos and tuple in file. 
  tree->close();}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analyseStepping(const G4Step* aStep)
{
  IHistogram1D*   histo_2 =dynamic_cast<IHistogram1D *> ( tree->find("2") );
  IHistogram1D*   histo_3 =dynamic_cast<IHistogram1D *> ( tree->find("3") );
  IHistogram1D*   histo_4 =dynamic_cast<IHistogram1D *> ( tree->find("4") );
  IHistogram1D*   histo_5 =dynamic_cast<IHistogram1D *> ( tree->find("5") );
  IHistogram1D*   histo_6 =dynamic_cast<IHistogram1D *> ( tree->find("6") );
  IHistogram1D*   histo_8 =dynamic_cast<IHistogram1D *> ( tree->find("8") );
  IHistogram1D*   histo_10 =dynamic_cast<IHistogram1D *> ( tree->find("10") );  
  G4double gammaAtTheDetPre=0;
  G4double protonsAtTheDetPre=0;
  G4double gammaLeavingSample=0;
  G4double gammaLeavingSampleSmeared=0;
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
	    /*gammaLeavingSampleSmeared = 
	      ResponseFunction(aStep->GetPreStepPoint()->GetKineticEnergy());
	      if(histo_10) {
	      histo_10->fill(gammaLeavingSampleSmeared/keV);*/
	    
	    gammaLeavingSample =(aStep->GetPreStepPoint()->GetKineticEnergy());
	    if(histo_4) {
	      histo_4->fill(gammaLeavingSample/keV);}
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
	      histo_2->fill(gammaBornInSample);
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
	      histo_3->fill(eleBornInSample);
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
		histo_6->fill( gammaAtTheDetPre);
	      }
	    }
	  else if ((aStep->GetTrack()->GetDynamicParticle()
		    ->GetDefinition()-> GetParticleName()) == "proton" ) 
	    {
	      protonsAtTheDetPre = 
		(aStep->GetPreStepPoint()->GetKineticEnergy());
	      if(histo_8) {
		histo_8->fill( protonsAtTheDetPre);
	      }
	    }
	}
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisManager::analyseEnergyDep(G4double energyDep)
{
 IHistogram1D* histo_1 =dynamic_cast<IHistogram1D *> ( tree->find("1") );
 
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
 IHistogram1D*   histo_7 =dynamic_cast<IHistogram1D *> ( tree->find("1") );
  histo_7->fill(energy/keV);
}





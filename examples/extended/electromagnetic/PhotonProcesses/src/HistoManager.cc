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
// $Id: HistoManager.cc,v 1.1 2004-04-28 11:12:39 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef USE_AIDA

#include "HistoManager.hh"
#include "HistoMessenger.hh"

#include "AIDA/AIDA.h"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:tree(0),hf(0),factoryOn(false)
{
  histoMessenger = new HistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ 
  if (factoryOn) { 
    tree->commit();       // Writing the histograms to the file
    tree->close();        // and closing the tree (and the file)
 
    delete hf; 
    delete tree;
    factoryOn = false;
  }

  delete histoMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetFactory(G4String fileName)
{ 
  if (factoryOn) {
    G4cout << "\n--->HistoManager::SetFactory(): factory already exists"
           << G4endl;
    return;
  }
    	   
 // Creating the analysis factory
 AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
  AIDA::ITreeFactory* tf = af->createTreeFactory();
 
 // Creating a tree mapped to an hbook file.
 G4bool readOnly  = false;
 G4bool createNew = true;
 tree = tf->create(fileName, "hbook", readOnly, createNew);

 // Creating a histogram factory, whose histograms will be handled by the tree
 hf = af->createHistogramFactory(*tree);
 
 // histograms
 for (G4int k=0; k<MaxHisto; k++) { histo[k] = 0; histoUnit[k] = 1.;}
  
 delete tf;
 delete af;
 factoryOn = true;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih, 
                 G4int nbBins, G4double valmin, G4double valmax, G4String unit)
{ 
  if (!factoryOn) SetFactory("photonprocesses.paw");
   
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }	 

  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6"};
  const G4String title[] = 
                { "dummy",						//0
                  "scattered primary particle: energy spectrum",	//1
                  "scattered primary particle: costheta distribution",	//2
                  "charged secondaries: energy spectrum",		//3
                  "charged secondaries: costheta distribution",		//4
                  "neutral secondaries: energy spectrum",		//5
                  "neutral secondaries: costheta distribution",		//6
                 };
		 
  G4String titl = title[ih];
  G4double vmin = valmin, vmax = valmax;
  
  if (unit != "none") {
    titl = title[ih] + " (" + unit + ")";
    histoUnit[ih] = G4UnitDefinition::GetValueOf(unit);
    vmin = valmin/histoUnit[ih]; vmax = valmax/histoUnit[ih];
  }
  
  histo[ih] = hf->createHistogram1D(id[ih],titl,nbBins,vmin,vmax);
  
  binWidth[ih] = (valmax-valmin)/nbBins;
  
  G4cout << "----> SetHisto " << ih << ": " << titl << ";  "
         << nbBins << " bins from " 
         << vmin << " " << unit << " to " << vmax << " " << unit << G4endl;   		 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

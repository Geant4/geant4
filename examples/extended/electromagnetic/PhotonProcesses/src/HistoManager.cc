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
// $Id: HistoManager.cc,v 1.2 2004/06/10 15:55:37 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
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
  fileName = "photonprocesses.paw";
  
  // histograms
  for (G4int k=0; k<MaxHisto; k++) { histo[k] = 0; exist[k] = false;}
   
  histoMessenger = new HistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ 
  SaveFactory();
  delete histoMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetFactory()
{ 
  if (factoryOn) {
    G4cout << "\n--->HistoManager::SetFactory(): factory already exists"
           << G4endl;
    SaveFactory();
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
 
 // create selected histograms
 for (G4int k=0; k<MaxHisto; k++) {
   if (exist[k]) histo[k] = hf->createHistogram1D( Label[k], Title[k],
                                                   Nbins[k], Vmin[k], Vmax[k]);
  }
        
 delete tf;
 delete af;
 factoryOn = true;  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SaveFactory()
{ 
  if (factoryOn) { 
    tree->commit();       // Writing the histograms to the file
    tree->close();        // and closing the tree (and the file)
 
    delete hf; 
    delete tree;
    factoryOn = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih, 
                 G4int nbins, G4double valmin, G4double valmax, G4String unit)
{   
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
                  "neutral secondaries: costheta distribution"		//6
                 };
		 
  G4String titl = title[ih];
  G4double vmin = valmin, vmax = valmax;
  Unit[ih] = 1.;
  
  if (unit != "none") {
    titl = title[ih] + " (" + unit + ")";
    Unit[ih] = G4UnitDefinition::GetValueOf(unit);
    vmin = valmin/Unit[ih]; vmax = valmax/Unit[ih];
  }
  
  exist[ih] = true;
  Label[ih] = id[ih];
  Title[ih] = titl;
  Nbins[ih] = nbins; 
  Vmin[ih]  = vmin; 
  Vmax[ih]  = vmax; 
  Width[ih] = (valmax-valmin)/nbins;
  
  G4cout << "----> SetHisto " << ih << ": " << titl << ";  "
         << nbins << " bins from " 
         << vmin << " " << unit << " to " << vmax << " " << unit << G4endl;
   		 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::RemoveHisto(G4int ih) 
{ 
 if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::RemoveHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }
  	  
  histo[ih] = 0; exist[ih] = false;     		 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

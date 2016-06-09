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
// $Id: HistoManager.cc,v 1.10 2005/06/01 10:48:40 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "HistoMessenger.hh"
#include "G4UnitsTable.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:af(0),tree(0),factoryOn(false)
{
#ifdef G4ANALYSIS_USE
  // Creating the analysis factory
  af = AIDA_createAnalysisFactory();
  if(!af) {
    G4cout << "TestEm1::HistoManager::HistoManager :" 
           << " problem creating the AIDA analysis factory."
           << G4endl;
  }	   
#endif
 
  fileName[0] = "testem3";
  fileType    = "hbook";  
  fileOption  = "--noErrors uncompress";  
  // histograms
  for (G4int k=0; k<MaxHisto; k++) { 
    histo[k] = 0;
    exist[k] = false;
    Unit[k]  = 1.0;
    Width[k] = 1.0;
  }   
  histoMessenger = new HistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ 
  delete histoMessenger;
  
#ifdef G4ANALYSIS_USE  
  delete af;
#endif     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
#ifdef G4ANALYSIS_USE
  if(!af) return;    	    
 
 // Creating a tree mapped to an hbook file.
 fileName[1] = fileName[0] + "." + fileType;
 G4bool readOnly  = false;
 G4bool createNew = true;
 AIDA::ITreeFactory* tf  = af->createTreeFactory(); 
 tree = tf->create(fileName[1], fileType, readOnly, createNew, fileOption);
 delete tf;
 if(!tree) {
   G4cout << "TestEm1::HistoManager::book :" 
          << " problem creating the AIDA tree with "
          << " storeName = " << fileName[1]
          << " storeType = " << fileType
          << " readOnly = "  << readOnly
          << " createNew = " << createNew
          << " options = "   << fileOption
          << G4endl;
   return;
 }
 
 // Creating a histogram factory, whose histograms will be handled by the tree
 AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);
 
 // create selected histograms
 for (G4int k=1; k<MaxHisto; k++) {
   if (exist[k]) {
     histo[k] = hf->createHistogram1D( Label[k], Title[k],
                                                 Nbins[k], Vmin[k], Vmax[k]);
     factoryOn = true;
   }  						  					   
 }
 delete hf;
 if (factoryOn) 
     G4cout << "\n----> Histogram Tree is opened in " << fileName[1]  << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
#ifdef G4ANALYSIS_USE
  if (factoryOn) {
    tree->commit();       // Writing the histograms to the file
    tree->close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved in " << fileName[1] << G4endl;

    delete tree;
    tree = 0;
    factoryOn = false;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
  if(exist[ih]) histo[ih]->fill(xbin/Unit[ih], weight);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih, 
            G4int nbins, G4double valmin, G4double valmax, const G4String& unit)
{   
  if (ih < 1 || ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }

 // histo 1 : energy deposit in absorber 1
 // histo 2 : energy deposit in absorber 2
 // ...etc...........
 // MaxAbsor = 10 (-1)
 // 
 // histo 11 : longitudinal profile of energy deposit in absorber 1 (MeV)
 // histo 12 : longitudinal profile of energy deposit in absorber 2 (MeV)  
 // ...etc...........  
 // 
 // histo 21 : forward  energy flow (MeV)
 // histo 22 : backward energy flow (MeV)  
  	 
  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
			 "20","21","22"};
			 
  G4String title;
  G4double vmin = valmin, vmax = valmax;
    			 
  if (ih < MaxAbsor) {			 
    title = "Edep in absorber " + id[ih] + " (" + unit + ")";
    Unit[ih] = G4UnitDefinition::GetValueOf(unit);   
    vmin = valmin/Unit[ih]; vmax = valmax/Unit[ih];  
  } else if (ih > MaxAbsor && ih < 2*MaxAbsor) {
    title = "longit. profile of Edep (MeV/event) in absorber " 
           + id[ih-MaxAbsor];
  } else if (ih == 2*MaxAbsor+1) {
    title = "Forward energy flow (MeV/event)";
  } else if (ih == 2*MaxAbsor+2) {
    title = "Backward energy flow (MeV/event)";
  } else return;
        
  exist[ih] = true;
  Label[ih] = id[ih];
  Title[ih] = title;
  Nbins[ih] = nbins; 
  Vmin[ih]  = vmin; 
  Vmax[ih]  = vmax; 
  Width[ih] = (valmax-valmin)/nbins;
  
  G4cout << "----> SetHisto " << ih << ": " << title << ";  "
         << nbins << " bins from " 
         << vmin << " " << unit << " to " << vmax << " " << unit << G4endl;
   		 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
  if(exist[ih]) histo[ih]->scale(fac);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::RemoveHisto(G4int ih) 
{ 
 if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::RemoveHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }
  	  
  histo[ih] = 0; exist[ih] = false;     		 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



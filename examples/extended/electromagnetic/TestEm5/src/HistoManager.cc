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
// $Id: HistoManager.cc,v 1.6 2004/02/19 18:18:53 maire Exp $
// GEANT4 tag $Name: geant4-06-01 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4ANALYSIS_USE

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
  if (!factoryOn) SetFactory("testem5.paw");
   
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }	 

  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16"};
  const G4String title[] = 
                { "dummy",						//0
                  "energy deposit in absorber",				//1
                  "energy of charged secondaries at creation",		//2
                  "energy of gammas at creation (log10(ekin/MeV))",	//3
		  "(transmit, charged) : kinetic energy at exit",	//4
		  "(transmit, charged) : space angle at exit",		//5
		  "(transmit, charged) : projected angle at exit",	//6
		  "(transmit, charged) : projected position at exit",	//7
		  "(transmit, neutral) : kinetic energy at exit",	//8
		  "(transmit, neutral) : space angle at exit",		//9
		  "(transmit, neutral) : projected angle at exit",	//10
		  "(reflect , charged) : kinetic energy at exit",	//11
		  "(reflect , charged) : space angle at exit",		//12
		  "(reflect , charged) : projected angle at exit",	//13
		  "(reflect , neutral) : kinetic energy at exit",	//14
		  "(reflect , neutral) : space angle at exit",		//15
		  "(reflect , neutral) : projected angle at exit",	//16
                 };
		 
  G4String titl = title[ih];
  G4double vmin = valmin, vmax = valmax;
  
  if (ih == 3) { vmin = log10(valmin/MeV); vmax = log10(valmax/MeV);}
  else if (unit != "none") {
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

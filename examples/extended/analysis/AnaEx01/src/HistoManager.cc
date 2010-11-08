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
// $Id: HistoManager.cc,v 1.1 2010-11-08 10:38:44 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:af(0),tree(0)
{
#ifdef G4ANALYSIS_USE
  // Creating the analysis factory
  //
  af = AIDA_createAnalysisFactory();
  if(!af) {
    G4cout << " HistoManager::HistoManager :" 
           << " problem creating the AIDA analysis factory."
           << G4endl;
  }	   
#endif
      
  // histograms
  for (G4int k=0; k<MaxHisto; k++) histo[k] = 0;
    
  // ntuple
  ntupl = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{  
#ifdef G4ANALYSIS_USE  
  delete af;
#endif     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
#ifdef G4ANALYSIS_USE
  if(!af) return;    	    
 
 // Creating a tree container to handle histograms and ntuples.
 // This tree is associated to an output file.
 //
 G4String fileName = "AnaEx01";
 G4String fileType    = "root";		// hbook  root  xml
 G4String fileOption  = " ";
 //// G4String fileOption  = "uncompress compress=no";		//for xml     
 //// G4String fileOption  = "--noErrors";			//for hbook

 fileName = fileName + "." + fileType;
 G4bool readOnly  = false;
 G4bool createNew = true;
 AIDA::ITreeFactory* tf  = af->createTreeFactory(); 
 tree = tf->create(fileName, fileType, readOnly, createNew, fileOption);
 delete tf;
 if(!tree) {
   G4cout << " HistoManager::book :" 
          << " problem creating the AIDA tree with "
          << " storeName = " << fileName
          << " storeType = " << fileType
          << " readOnly = "  << readOnly
          << " createNew = " << createNew
          << " options = "   << fileOption
          << G4endl;
   return;
 }
 
 // Creating a histogram factory, whose histograms will be handled by the tree
 //
 AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);
 
 // create histos in subdirectory "histograms"
 //
 tree->mkdir("histograms");
 tree->cd("histograms");
  
 histo[1] = hf->createHistogram1D("1", "Edep in absorber", 100, 0., 800*MeV);
 if (!histo[1]) G4cout << "\n can't create histo 1" << G4endl;
 histo[2] = hf->createHistogram1D("2", "Edep in gap", 100, 0., 100*MeV);
 if (!histo[2]) G4cout << "\n can't create histo 2" << G4endl;
 histo[3] = hf->createHistogram1D("3", "trackL in absorber", 100, 0., 1*m);
 if (!histo[3]) G4cout << "\n can't create histo 3" << G4endl;
 histo[4] = hf->createHistogram1D("4", "trackL in gap", 100, 0., 50*cm);
 if (!histo[4]) G4cout << "\n can't create histo 4" << G4endl;  

 delete hf;
 tree->cd(".."); 
 
 // Creating a ntuple factory, handled by the tree
 //
 AIDA::ITupleFactory* ntf = af->createTupleFactory(*tree);
 
 // create 1 ntuple in subdirectory "tuples"
 //
 tree->mkdir("tuples");
 tree->cd("tuples");
  
 ntupl = ntf->create("101", "Edep and TrackL", "double Eabs, Egap, Labs, Lgap");
 
 delete ntf;
 tree->cd("..");
    
 G4cout << "\n----> Histogram Tree is opened in " << fileName << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
#ifdef G4ANALYSIS_USE
  if (af && tree) {
    tree->commit();       // Writing the histograms to the file
    tree->close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;

    delete tree;
    tree = 0;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << " does not exist. (xbin=" << xbin << " weight=" << weight << ")"
	   << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
 if (histo[ih]) histo[ih]->fill(xbin, weight);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << " does not exist. (fac=" << fac << ")" << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
   if (histo[ih]) histo[ih]->scale(fac);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4int column, G4double value)
{
  if (column > 3) {
    G4cout << "---> warning from HistoManager::FillNtuple : " 
     << "column=" << column << " value=" << value << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
  if (ntupl) ntupl->fill(column, value);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddRowNtuple()
{
#ifdef G4ANALYSIS_USE
  if (ntupl) ntupl->addRow();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
#ifdef G4ANALYSIS_USE
  if(histo[1]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
       << " EAbs : mean = " << G4BestUnit(histo[1]->mean(), "Energy") 
               << " rms = " << G4BestUnit(histo[1]->rms(),  "Energy") << G4endl;
    G4cout 	       
       << " EGap : mean = " << G4BestUnit(histo[2]->mean(), "Energy") 
               << " rms = " << G4BestUnit(histo[2]->rms(),  "Energy") << G4endl;
    G4cout 
       << " LAbs : mean = " << G4BestUnit(histo[3]->mean(), "Length") 
               << " rms = " << G4BestUnit(histo[3]->rms(),  "Length") << G4endl;
    G4cout 
       << " LGap : mean = " << G4BestUnit(histo[4]->mean(), "Length") 
               << " rms = " << G4BestUnit(histo[4]->rms(),  "Length") << G4endl;

  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



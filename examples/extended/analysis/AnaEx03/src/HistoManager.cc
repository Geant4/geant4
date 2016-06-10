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
/// \file analysis/AnaEx03/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 74272 2013-10-02 14:48:50Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:fAF(0),fTree(0), fNtuple1(0), fNtuple2(0)
{
#ifdef G4ANALYSIS_USE
  // Creating the analysis factory
  //
  fAF = AIDA_createAnalysisFactory();
  if(!fAF) {
    G4cout << " HistoManager::HistoManager :" 
           << " problem creating the AIDA analysis factory."
           << G4endl;
  }           
#endif
      
  // histograms
  for (G4int k=0; k<MaxHisto; k++) fHisto[k] = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{  
#ifdef G4ANALYSIS_USE  
  delete fAF;
#endif     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
#ifdef G4ANALYSIS_USE
  if(!fAF) return;                
 
 // Creating a tree container to handle histograms and ntuples.
 // This tree is associated to an output file.
 //
 G4String fileName = "AnaEx03";
 G4String fileType    = "root";                // hbook  root  xml
 G4String fileOption  = " ";
 //// G4String fileOption  = "uncompress compress=no";  //for xml     
 //// G4String fileOption  = "--noErrors";              //for hbook

 fileName = fileName + "." + fileType;
 G4bool readOnly  = false;
 G4bool createNew = true;
 AIDA::ITreeFactory* tf  = fAF->createTreeFactory(); 
 fTree = tf->create(fileName, fileType, readOnly, createNew, fileOption);
 delete tf;
 if(!fTree) {
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
 AIDA::IHistogramFactory* hf = fAF->createHistogramFactory(*fTree);
 
 // create histos in subdirectory "histograms"
 //
 fTree->mkdir("histograms");
 fTree->cd("histograms");
  
 fHisto[1] = hf->createHistogram1D("1", "Edep in absorber", 100, 0., 800*MeV);
 if (!fHisto[1]) G4cout << "\n can't create histo 1" << G4endl;
 fHisto[2] = hf->createHistogram1D("2", "Edep in gap", 100, 0., 100*MeV);
 if (!fHisto[2]) G4cout << "\n can't create histo 2" << G4endl;
 fHisto[3] = hf->createHistogram1D("3", "trackL in absorber", 100, 0., 1*m);
 if (!fHisto[3]) G4cout << "\n can't create histo 3" << G4endl;
 fHisto[4] = hf->createHistogram1D("4", "trackL in gap", 100, 0., 50*cm);
 if (!fHisto[4]) G4cout << "\n can't create histo 4" << G4endl;  

 delete hf;
 fTree->cd(".."); 
 
 // Creating a ntuple factory, handled by the tree
 //
 AIDA::ITupleFactory* ntf = fAF->createTupleFactory(*fTree);
 
 // create 1 ntuple in subdirectory "tuples"
 //
 fTree->mkdir("tuples");
 fTree->cd("tuples");
  
 fNtuple1 = ntf->create("101", "Edep", "double Eabs, Egap");
 fNtuple2 = ntf->create("102", "TrackL", "double Labs, Lgap");
 
 delete ntf;
 fTree->cd("..");
    
 G4cout << "\n----> Histogram Tree is opened in " << fileName << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
#ifdef G4ANALYSIS_USE
  if (fAF && fTree) {
    fTree->commit();       // Writing the histograms to the file
    fTree->close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;

    delete fTree;
    fTree = 0;
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
 if (fHisto[ih]) fHisto[ih]->fill(xbin, weight);
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
   if (fHisto[ih]) fHisto[ih]->scale(fac);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4ANALYSIS_USE
void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap,
                              G4double trackLAbs, G4double trackLGap)
{                              
  if (fNtuple1) {
    fNtuple1->fill(0, energyAbs);
    fNtuple1->fill(1, energyGap);
    fNtuple1->addRow();
  }  
  if (fNtuple2) {
    fNtuple2->fill(0, trackLAbs);
    fNtuple2->fill(1, trackLGap);
    fNtuple2->addRow();
  }  
}
#else
void HistoManager::FillNtuple(G4double, G4double, G4double, G4double)
{
  return;
}
#endif

void HistoManager::PrintStatistic()
{
#ifdef G4ANALYSIS_USE
  if(fHisto[1]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
      << " EAbs : mean = " << G4BestUnit(fHisto[1]->mean(), "Energy") 
              << " rms = " << G4BestUnit(fHisto[1]->rms(),  "Energy") << G4endl;
    G4cout                
      << " EGap : mean = " << G4BestUnit(fHisto[2]->mean(), "Energy") 
              << " rms = " << G4BestUnit(fHisto[2]->rms(),  "Energy") << G4endl;
    G4cout 
      << " LAbs : mean = " << G4BestUnit(fHisto[3]->mean(), "Length") 
              << " rms = " << G4BestUnit(fHisto[3]->rms(),  "Length") << G4endl;
    G4cout 
      << " LGap : mean = " << G4BestUnit(fHisto[4]->mean(), "Length") 
              << " rms = " << G4BestUnit(fHisto[4]->rms(),  "Length") << G4endl;

  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



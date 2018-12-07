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
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "AIDA/AIDA.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:fAF(0),fTree(0), fNtuple1(0), fNtuple2(0)
{
  // Creating the analysis factory
  //
  fAF = AIDA_createAnalysisFactory();
  if(!fAF) {
    G4cout << " HistoManager::HistoManager :" 
           << " problem creating the AIDA analysis factory."
           << G4endl;
  }           
      
  // histograms
  for (G4int k=0; k<kMaxHisto; k++) fHisto[k] = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{  
  delete fAF;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{ 
  if (! fAF) return;                
 
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
  
  // id = 0
  fHisto[0] = hf->createHistogram1D("EAbs", "EAbs: Edep in absorber", 100, 0., 800*MeV);
  // id = 1
  fHisto[1] = hf->createHistogram1D("EGap", "EGap: Edep in gap", 100, 0., 100*MeV);
  // id = 2
  fHisto[2] = hf->createHistogram1D("LAbs", "LAbs: trackL in absorber", 100, 0., 1*m);
  // id = 3
  fHisto[3] = hf->createHistogram1D("LGap", "LGap: trackL in gap", 100, 0., 50*cm);

  for ( G4int i=0; i<kMaxHisto; ++i ) {
    if (! fHisto[i]) G4cout << "\n can't create histo " << i << G4endl;
  }  

  delete hf;
  fTree->cd(".."); 
 
  // Creating a ntuple factory, handled by the tree
  //
  AIDA::ITupleFactory* ntf = fAF->createTupleFactory(*fTree);
 
  // create 1 ntuple in subdirectory "tuples"
  //
  fTree->mkdir("tuples");
  fTree->cd("tuples");
  
  fNtuple1 = ntf->create("Ntuple1", "Edep", "double Eabs, Egap");
  fNtuple2 = ntf->create("Ntuple2", "TrackL", "double Labs, Lgap");
 
  delete ntf;
  fTree->cd("..");
    
  G4cout << "\n----> Output file is open in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save()
{ 
  if (! (fAF && fTree)) return;

  fTree->commit();       // Writing the histograms to the file
  fTree->close();        // and closing the tree (and the file)
  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;

  delete fTree;
  fTree = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= kMaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << " does not exist. (xbin=" << xbin << " weight=" << weight << ")"
           << G4endl;
    return;
  }

  if (fHisto[ih]) fHisto[ih]->fill(xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= kMaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << " does not exist. (fac=" << fac << ")" << G4endl;
    return;
  }

  if (fHisto[ih]) fHisto[ih]->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

void HistoManager::PrintStatistic()
{
  G4cout << "\n ----> print histograms statistic \n" << G4endl;
  for ( G4int i=0; i<kMaxHisto; ++i ) {
    AIDA::IHistogram1D* h1 = fHisto[i];
    G4String title = h1->title();  
    // extract name as first 4 characters from title, as aida seems not to keep
    // histogram name
    const G4String name = title(0,4);  

    G4String unitCategory;
    if (name[0] == 'E' ) unitCategory = "Energy"; 
    if (name[0] == 'L' ) unitCategory = "Length";

    G4cout << name
           << ": mean = " << G4BestUnit(h1->mean(), unitCategory) 
           << " rms = " << G4BestUnit(h1->rms(), unitCategory ) 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

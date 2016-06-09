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
// $Id: HistoManager.cc 48195 2011-02-02 15:33:39Z jjacquem $
// GEANT4 tag $Name: geant4-09-04 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

#ifdef G4ANALYSIS_USE
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:rootFile(0),ntupl(0), Eabs(0), Egap(0) ,Labs(0), Lgap(0)
{
      
  // histograms
  for (G4int k=0; k<MaxHisto; k++) histo[k] = 0;
    
  // ntuple
  ntupl = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
#ifdef G4ANALYSIS_USE  
    if ( rootFile ) delete rootFile;
#endif    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
#ifdef G4ANALYSIS_USE
 
 // Creating a tree container to handle histograms and ntuples.
 // This tree is associated to an output file.
 //
 G4String fileName = "AnaEx02.root";
 rootFile = new TFile(fileName,"RECREATE");
 if(!rootFile) {
   G4cout << " HistoManager::book :" 
          << " problem creating the ROOT TFile "
          << G4endl;
   return;
 }
   
 histo[1] = new TH1D("1", "Edep in absorber", 100, 0., 800*MeV);
 if (!histo[1]) G4cout << "\n can't create histo 1" << G4endl;
 histo[2] = new TH1D("2", "Edep in gap", 100, 0., 100*MeV);
 if (!histo[2]) G4cout << "\n can't create histo 2" << G4endl;
 histo[3] = new TH1D("3", "trackL in absorber", 100, 0., 1*m);
 if (!histo[3]) G4cout << "\n can't create histo 3" << G4endl;
 histo[4] = new TH1D("4", "trackL in gap", 100, 0., 50*cm);
 if (!histo[4]) G4cout << "\n can't create histo 4" << G4endl;  

 // create 1 ntuple in subdirectory "tuples"
 //
 ntupl = new TTree("101", "Edep and TrackL");
 ntupl->Branch("Eabs", &Eabs, "Eabs/D");
 ntupl->Branch("Egap", &Egap, "Egap/D");
 ntupl->Branch("Labs", &Labs, "Labs/D");
 ntupl->Branch("Lgap", &Lgap, "Lgap/D");

 
 G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
#ifdef G4ANALYSIS_USE
  if (rootFile) {
    rootFile->Write();       // Writing the histograms to the file
    rootFile->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;
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
 if  (histo[ih]) { histo[ih]->Fill(xbin, weight); }
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
   if (histo[ih]) histo[ih]->Scale(fac);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4double EnergyAbs, G4double EnergyGap,
                              G4double TrackLAbs , G4double TrackLGap )
{
 Eabs = EnergyAbs;
 Egap = EnergyGap;
 Labs = TrackLAbs;
 Lgap = TrackLGap;

#ifdef G4ANALYSIS_USE
  if (ntupl) ntupl->Fill();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
#ifdef G4ANALYSIS_USE
  if(histo[1]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " EAbs : mean = " << G4BestUnit(histo[1]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[1]->GetRMS(),  "Energy") << G4endl;
    G4cout 	       
    << " EGap : mean = " << G4BestUnit(histo[2]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[2]->GetRMS(),  "Energy") << G4endl;
    G4cout 
    << " LAbs : mean = " << G4BestUnit(histo[3]->GetMean(), "Length") 
            << " rms = " << G4BestUnit(histo[3]->GetRMS(),  "Length") << G4endl;
    G4cout 
    << " LGap : mean = " << G4BestUnit(histo[4]->GetMean(), "Length") 
            << " rms = " << G4BestUnit(histo[4]->GetRMS(),  "Length") << G4endl;

  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



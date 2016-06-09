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
/// \file analysis/AnaEx02/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 65827 2012-11-29 10:21:47Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

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
  if ( rootFile ) delete rootFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
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
   
 histo[1] = new TH1D("1", "Edep in absorber", 100, 0., 800*CLHEP::MeV);
 if (!histo[1]) G4cout << "\n can't create histo 1" << G4endl;
 histo[2] = new TH1D("2", "Edep in gap", 100, 0., 100*CLHEP::MeV);
 if (!histo[2]) G4cout << "\n can't create histo 2" << G4endl;
 histo[3] = new TH1D("3", "trackL in absorber", 100, 0., 1*CLHEP::m);
 if (!histo[3]) G4cout << "\n can't create histo 3" << G4endl;
 histo[4] = new TH1D("4", "trackL in gap", 100, 0., 50*CLHEP::cm);
 if (!histo[4]) G4cout << "\n can't create histo 4" << G4endl;  

 // create 1 ntuple in subdirectory "tuples"
 //
 ntupl = new TTree("101", "Edep and TrackL");
 ntupl->Branch("Eabs", &Eabs, "Eabs/D");
 ntupl->Branch("Egap", &Egap, "Egap/D");
 ntupl->Branch("Labs", &Labs, "Labs/D");
 ntupl->Branch("Lgap", &Lgap, "Lgap/D");

 
 G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
  if (rootFile) {
    rootFile->Write();       // Writing the histograms to the file
    rootFile->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;
  }
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
 if  (histo[ih]) { histo[ih]->Fill(xbin, weight); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << " does not exist. (fac=" << fac << ")" << G4endl;
    return;
  }
  if (histo[ih]) histo[ih]->Scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap,
                              G4double trackLAbs , G4double trackLGap )
{
 Eabs = energyAbs;
 Egap = energyGap;
 Labs = trackLAbs;
 Lgap = trackLGap;

  if (ntupl) ntupl->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



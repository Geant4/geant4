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
/// \file electromagnetic/TestEm9/src/Histo.cc
/// \brief Implementation of the Histo class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   Histo - Generic histogram/ntuple manager class
//
//
// Author:      V.Ivanchenko 30.10.03
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Histo.hh"
#include "HistoMessenger.hh"
#include "G4RootAnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Histo::Histo()
 : fManager(0),
   fMessenger(0)
{
  fMessenger = new HistoMessenger(this);
 
  fHistName   = "test";
  fHistType   = "root";
  fTupleName  = "tuple";
  fTupleTitle = "test";
  fNHisto     = 0;
  fVerbose    = 0;
  fDefaultAct = true;
  fHistoActive= false;
  fNtupleActive= false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Histo::~Histo()
{
  delete fMessenger;
  delete fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::Book()
{
  if(!(fHistoActive || fNtupleActive)) { return; }

  // Always creating analysis manager
  fManager = G4RootAnalysisManager::Instance(); 

  // Creating a tree mapped to a new hbook file.
  G4String nam = fHistName + "." + fHistType;

  // Open file histogram file
  if(!fManager->OpenFile(nam)) {
    G4cout << "Histo::Book: ERROR open file <" << nam << ">" << G4endl;
    fHistoActive = false;
    fNtupleActive = false;
    return; 
  }
  G4cout << "### Histo::Save: Opended file <" << nam << ">  for " 
         << fNHisto << " histograms " << G4endl;

  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<fNHisto; ++i) {
    if(fActive[i]) {
      G4String ss = "h" + fIds[i];
      fHisto[i] = 
       fManager->CreateH1(ss, fTitles[i], fBins[i], fXmin[i], fXmax[i]);
      if(fVerbose > 0) {
        G4cout << "Created histogram #" << i << "  id= " << fHisto[i]
               << "  "  << ss << "  " << fTitles[i] << G4endl;
      }
    }
  }
  // Creating a tuple factory, whose tuples will be handled by the tree
  if(fNtupleActive) {
    fManager->CreateNtuple(fTupleName,fTupleTitle); 
    G4int i;
    G4int n = fNtupleI.size();
    for(i=0; i<n; ++i) { 
      if(fTupleI[i] == -1) {  
       fTupleI[i] = fManager->CreateNtupleIColumn(fNtupleI[i]); 
      }
    }
    n = fNtupleF.size();
    for(i=0; i<n; ++i) { 
      if(fTupleF[i] == -1) {  
       fTupleF[i] = fManager->CreateNtupleFColumn(fNtupleF[i]); 
      }
    }
    n = fNtupleD.size();
    for(i=0; i<n; ++i) { 
      if(fTupleD[i] == -1) {  
       fTupleD[i] = fManager->CreateNtupleDColumn(fNtupleD[i]); 
      }
    }
  }
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::Save()
{
  if(!(fHistoActive || fNtupleActive)) { return; }

  // Creating a tree mapped to a new hbook file.
  G4String nam = fHistName + "." + fHistType;

  // Write histogram file
  if(!fManager->Write()) {
    G4Exception ("Histo::Save()", "hist01", FatalException, 
                 "Cannot write ROOT file.");
  }
  if(fVerbose > 0) {
    G4cout << "### Histo::Save: Histograms and Ntuples are saved" << G4endl;
  }
  if(fManager->CloseFile() && fVerbose > 0) {
    G4cout << "                 File is closed" << G4endl;
  }
  delete G4RootAnalysisManager::Instance();
  fManager = 0;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::Add1D(const G4String& id, const G4String& name, G4int nb, 
                  G4double x1, G4double x2, G4double u)
{
  if(fVerbose > 0) {
    G4cout << "Histo::Add1D: New histogram will be booked: #" << id 
           << "  <" << name 
           << "  " << nb << "  " << x1 << "  " << x2 << "  " << u 
           << G4endl;
  }
  ++fNHisto;
  x1 /= u;
  x2 /= u;
  fActive.push_back(fDefaultAct);
  fBins.push_back(nb);
  fXmin.push_back(x1);
  fXmax.push_back(x2);
  fUnit.push_back(u);
  fIds.push_back(id);
  fTitles.push_back(name);
  fHisto.push_back(-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::SetHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double u)
{
  if(i>=0 && i<fNHisto) {
    if(fVerbose > 0) {
      G4cout << "Histo::SetHisto1D: #" << i  
             << "  " << nb << "  " << x1 << "  " << x2 << "  " << u 
             << G4endl;
    }
    fBins[i] = nb;
    fXmin[i] = x1;
    fXmax[i] = x2;
    fUnit[i] = u;
    fActive[i] = true;
    fHistoActive = true;
  } else {
    G4cout << "Histo::SetHisto1D: WARNING! wrong histogram index " 
          << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::Activate(G4int i, G4bool val)
{
  if(fVerbose > 1) {
    G4cout << "Histo::Activate: Histogram: #" << i << "   "  
           << val << G4endl;   
  }
  if(i>=0 && i<fNHisto) { 
    fActive[i] = val; 
    if(val) { fHistoActive = true; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::Fill(G4int i, G4double x, G4double w)
{
  if(!fHistoActive) { return; }
  if(fVerbose > 1) {
    G4cout << "Histo::Fill: Histogram: #" << i << " at x= " << x 
           << "  weight= " << w
           << G4endl;   
  }
  if(i>=0 && i<fNHisto) {
    if(fActive[i]) { fManager->FillH1(fHisto[i], x/fUnit[i], w); }
  } else {
    G4cout << "Histo::Fill: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::ScaleH1(G4int i, G4double x)
{
  if(!fHistoActive) { return; }
  if(fVerbose > 0) {
    G4cout << "Histo::Scale: Histogram: #" 
          << i << " by factor " << x << G4endl;   
  }
  if(i>=0 && i<fNHisto) {
    if(fActive[i]) { fManager->GetH1(fHisto[i])->scale(x); }
  } else {
    G4cout << "Histo::Scale: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::AddTuple(const G4String& w1)
{
  fTupleTitle = w1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::AddTupleI(const G4String& w1)
{
  fNtupleActive = true;
  fNtupleI.push_back(w1);
  fTupleI.push_back(-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::AddTupleF(const G4String& w1)
{
  fNtupleActive = true;
  fNtupleF.push_back(w1);
  fTupleF.push_back(-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::AddTupleD(const G4String& w1)
{
  fNtupleActive = true;
  fNtupleD.push_back(w1);
  fTupleD.push_back(-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::FillTupleI(G4int i, G4int x)
{
  if(!fNtupleActive) { return; }
  G4int n = fNtupleI.size();
  if(i >= 0 && i < n) {
    if(fVerbose > 1) {
      G4cout << "Histo::FillTupleI: i= " << i << "  id= " << fTupleI[i]
             << "   <" << fNtupleI[i] << "> = " << x << G4endl; 
    }
    fManager->FillNtupleIColumn(fTupleI[i], x); 
  } else {
    G4cout << "Histo::FillTupleI: WARNING! wrong ntuple index " 
           << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::FillTupleF(G4int i, G4float x)
{
  if(!fNtupleActive) { return; }
  G4int n = fNtupleF.size();
  if(i >= 0 && i < n) {
    if(fVerbose > 1) {
      G4cout << "Histo::FillTupleF: i= " << i << "  id= " << fTupleF[i]
             << "   <" << fNtupleF[i] << "> = " << x << G4endl; 
    }
    fManager->FillNtupleFColumn(fTupleF[i], x); 
  } else {
    G4cout << "Histo::FillTupleF: WARNING! wrong ntuple index " 
           << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::FillTupleD(G4int i, G4double x)
{
  if(!fNtupleActive) { return; }
  G4int n = fNtupleD.size();
  if(i >= 0 && i < n) {
    if(fVerbose > 1) {
      G4cout << "Histo::FillTupleD: i= " << i << "  id= " << fTupleD[i]
             << "   <" << fNtupleD[i] << "> = " << x << G4endl; 
    }
    fManager->FillNtupleDColumn(fTupleD[i], x); 
  } else {
    G4cout << "Histo::FillTupleD: WARNING! wrong ntuple index " 
           << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::AddRow()
{
  if(!fNtupleActive) { return; }
  fManager->AddNtupleRow();
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::SetFileName(const G4String& nam) 
{
  fHistName = nam;
  fHistoActive = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Histo::SetFileType(const G4String& nam) 
{
  // format other than ROOT is not tested
  if(nam == "root" || nam == "ROOT" )   { fHistType = "root"; }
  else if(nam == "xml" || nam == "XML") { fHistType = "xml"; }
  else if(nam == "ascii" || nam == "ASCII" || 
          nam == "Csv" || nam == "csv" || nam == "CSV") 
    { fHistType = "ascii"; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


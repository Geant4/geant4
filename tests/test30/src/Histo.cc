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
//---------------------------------------------------------------------------
//
// ClassName:   Histo - Generic histogram/ntuple manager class
//
//
// Author:      V.Ivanchenko 30.10.03
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Histo.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

#ifdef G4ANALYSIS_USE_ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::Histo()
{
  verbose    = 0;
  histName   = "test30";
  histType   = "hbook";
  nHisto     = 0;
  defaultAct = true;
  tupleName  = "tuple";
  tupleId    = "100";
  tupleList  = "";
  ntup       = 0;
  m_ROOT_file= 0;

#ifdef G4ANALYSIS_USE
  tree = 0;
  af   = 0; 
#endif
  //
#ifdef G4ANALYSIS_USE_ROOT
  histType   = "root";
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::~Histo()
{
#ifdef G4ANALYSIS_USE
  delete af;
#endif
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::book()
{
  G4String nam = histName + "." + histType;
  G4cout << "### Histo books " << nHisto 
	 << " histograms in <" << nam << ">" << G4endl;

#ifdef G4ANALYSIS_USE
  // Creating the analysis factory
  if(!af) af = AIDA_createAnalysisFactory();
  // Creating the tree factory
  AIDA::ITreeFactory* tf = af->createTreeFactory();

  // Creating a tree mapped to a new hbook file.
  tree = tf->create(nam,histType,false,true,"--noErrors uncompress");
  delete tf;
  if(tree) {
    G4cout << "Tree store  : <" << tree->storeName() << ">" << G4endl;
  } else {
    G4cout << "ERROR: Tree store " << histName  << " is not created!" << G4endl;
    return;
  }
  // Creating a histogram factory, whose histograms will be handled by the tree
  AIDA::IHistogramFactory* hf = af->createHistogramFactory( *tree );

  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<nHisto; i++) {
    if(active[i]) {
      histo[i] = hf->createHistogram1D(ids[i], titles[i], bins[i], xmin[i], xmax[i]);
    }
  }
  delete hf;
  // Creating a tuple factory, whose tuples will be handled by the tree
  if(tupleList != "") {
     AIDA::ITupleFactory* tpf = af->createTupleFactory( *tree );
     ntup = tpf->create(tupleId, tupleName, tupleList);
     delete tpf;
  }
#endif

  // Creating the ROOT file, trees, histos
#ifdef G4ANALYSIS_USE_ROOT
  m_ROOT_file = 
    new TFile(nam,"RECREATE","ROOT file with trees and histograms");
  if(m_ROOT_file)
    G4cout << "[Histo::book] File created: " << nam << G4endl;
  else
    G4Exception("[Histo::book] ERROR: file " + nam + " has not been created!");


  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<nHisto; i++) {
    if(active[i]) {
      G4String r_name = "h" + ids[i];
      m_ROOT_histo[i] = new TH1D(r_name, titles[i], bins[i], xmin[i], xmax[i]);
    } 
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::save()
{
#ifdef G4ANALYSIS_USE
  // Write histogram file
  if(tree) {
    tree->commit();
    G4cout << "Closing the tree..." << G4endl;
    tree->close();
    G4cout << "Histograms and Ntuples are saved" << G4endl;
    delete tree;
    tree = 0;
  }
#endif

  // Writing and closing the ROOT file
#ifdef G4ANALYSIS_USE_ROOT
  G4String nam = histName + "." + histType;
  G4cout << "[Histo::save] ROOT: file writing <" <<nam << ">"<< G4endl;
  m_ROOT_file->Write();
  G4cout << "[Histo::save] ROOT: files closing..." << G4endl;
  m_ROOT_file->Close();
  delete m_ROOT_file;
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::add1D(const G4String& id, const G4String& name, G4int nb, 
                  G4double x1, G4double x2, G4double u)
{
  if(verbose > 0) {
    G4cout << "New histogram will be booked: #" << id << "  <" << name 
           << "  " << nb << "  " << x1 << "  " << x2 << "  " << u 
           << G4endl;
  }
  nHisto++;
  x1 /= u;
  x2 /= u;
  active.push_back(defaultAct);
  bins.push_back(nb);
  xmin.push_back(x1);
  xmax.push_back(x2);
  unit.push_back(u);
  ids.push_back(id);
  titles.push_back(name);
  histo.push_back(0);
  m_ROOT_histo.push_back(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double u)
{
  if(i>=0 && i<nHisto) {
    if(verbose > 0) {
      G4cout << "Update histogram: #" << i  
             << "  " << nb << "  " << x1 << "  " << x2 << "  " << u 
             << G4endl;
    }
    bins[i] = nb;
    xmin[i] = x1;
    xmax[i] = x2;
    unit[i] = u;
  } else {
    G4cout << "Histo::setHisto1D: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::fill(G4int i, G4double x, G4double w)
{
  if(verbose > 1) {
    G4cout << "fill histogram: #" << i << " at x= " << x 
           << "  weight= " << w
           << G4endl;   
  }
  if(i>=0 && i<nHisto) {
    if(active[i]) {
#ifdef G4ANALYSIS_USE  
      if(histo[i]) histo[i]->fill((x/unit[i]), w);
#endif
      //
#ifdef G4ANALYSIS_USE_ROOT
      if(m_ROOT_histo[i]) m_ROOT_histo[i]->Fill(x/unit[i], w);
#endif
    }
  } else {
    if(verbose > 0) 
      G4cout << "Histo::fill: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::scale(G4int i, G4double x)
{
  if(verbose > 0) {
    G4cout << "Scale histogram: #" << i << " by factor " << x << G4endl;   
  }
  if(i>=0 && i<nHisto) {
    if(active[i]) {
#ifdef G4ANALYSIS_USE  
      if(histo[i]) histo[i]->scale(x);
#endif
      //
#ifdef G4ANALYSIS_USE_ROOT  
      if(m_ROOT_histo[i]) m_ROOT_histo[i]->Scale(x);
#endif
    }
  } else {
    G4cout << "Histo::scale: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::activate(G4int i, G4bool val)
{
  if(i>=0 && i<nHisto) {
    active[i] = val;
  } else {
    if(verbose > 0) 
      G4cout << "Histo::activate: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::addTuple(const G4String& w1, const G4String& w2, const G4String& w3)
{
  tupleId = w1;
  tupleName = w2;
  tupleList = w3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::fillTuple(const G4String& parname, G4double x)
{
  if(verbose > 1) {
    G4cout << "fill tuple by parameter <" << parname << "> = " << x << G4endl; 
  }
#ifdef G4ANALYSIS_USE  
  if(ntup) ntup->fill(ntup->findColumn(parname), (float)x);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::addRow()
{
#ifdef G4ANALYSIS_USE
  if(ntup) ntup->addRow();
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setFileName(const G4String& nam) 
{
  histName = nam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setFileType(const G4String& nam) 
{
  if(nam == "root" || nam == "hbook" || nam == "aida") histType = nam;
  else if(nam == "xml" || nam == "XML") histType = "aida";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


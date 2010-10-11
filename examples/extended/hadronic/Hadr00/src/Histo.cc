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
#include "HistoMessenger.hh"

#ifdef G4ANALYSIS_USE
#include <AIDA/AIDA.h>
#endif

Histo::Histo(G4int ver)
{
  verbose    = ver;
  histName   = "histo";
  histType   = "root";
  option     = "";
  nHisto     = 0;
  defaultAct = 1;
  messenger  = new HistoMessenger(this);
#ifdef G4ANALYSIS_USE
  tree = 0;
  af   = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::~Histo()
{
  delete messenger;
#ifdef G4ANALYSIS_USE
  delete af;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::book()
{
#ifdef G4ANALYSIS_USE
  G4cout << "### Histo books " << nHisto << " histograms " << G4endl;
  // Creating the analysis factory
  AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
  if(verbose>0) {
    G4cout<<"Histo books analysis factory"<<G4endl;
  }
  // Creating the tree factory
  AIDA::ITreeFactory* tf = af->createTreeFactory();
  if(verbose>0)
    G4cout<<"Histo books tree factory ......... "<<G4endl;

  // Creating a tree mapped to a new hbook file.
  G4String name = histName + "." + histType;
  tree = tf->create(name, histType, false, true, option);
  G4cout << "Histo: tree store : " << tree->storeName() << G4endl;

  // Creating a histogram factory, whose histograms will be handled by the tree
  AIDA::IHistogramFactory* hf = af->createHistogramFactory( *tree );

  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<nHisto; i++) {
    if(active[i]) {
      if(verbose>0) {
	G4cout<<"Histo::book: histogram "<< i << " id= " << ids[i] <<G4endl;
      }
      histo[i] = hf->createHistogram1D(ids[i], tittles[i], bins[i], xmin[i], xmax[i]);
    }
  }
  delete hf;
  delete tf;
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::save()
{
#ifdef G4ANALYSIS_USE
  // Write histogram file
  tree->commit();
  G4cout << "Closing the tree..." << G4endl;
  tree->close();
  G4cout << "Histograms and Ntuples are saved" << G4endl;
  delete tree;
  tree = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setFileType(const G4String& nam) 
{
  if(nam == "hbook" || nam == "root" || nam == "aida") { histType = nam; }
  else if(nam == "XML" || nam == "xml") { histType = "aida"; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::add1D(const G4String& id, const G4String& name, G4int nb,
                  G4double x1, G4double x2, G4double u)
{
  if(nHisto > 0) {
    for(G4int i=0; i<nHisto; ++i) {
      if(ids[i] == id) return;
    }
  }

  if(verbose > 0) {
    G4cout << "New histogram will be booked: #" << id << "  <" << name
           << "  " << nb << "  " << x1 << "  " << x2 << "  " << u
           << G4endl;
  }
  ++nHisto;
  x1 /= u;
  x2 /= u;
  active.push_back(defaultAct);
  bins.push_back(nb);
  xmin.push_back(x1);
  xmax.push_back(x2);
  unit.push_back(u);
  ids.push_back(id);
  tittles.push_back(name);
  histo.push_back(0);
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
    G4cout << "nHisto = " << nHisto << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::fill(G4int i, G4double x, G4double w)
{
  if(verbose > 1) {
    G4cout << "fill histogram: #" << i << " at x= " << x 
           << "  weight= " << w << " unit= " << unit[i]
           << G4endl;
  }   
#ifdef G4ANALYSIS_USE  
  if(i>=0 && i<nHisto) {
    histo[i]->fill(x/unit[i], w);
    //histo[i]->fill((float)(x/unit[i]), (float)w);
  } else {
    G4cout << "Histo::fill: WARNING! wrong histogram index " << i << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::scale(G4int i, G4double x)
{
  if(verbose > 0) {
    G4cout << "Scale histogram: #" << i << " by factor " << x << G4endl;   
  }
#ifdef G4ANALYSIS_USE  
  if(i>=0 && i<nHisto) {
    histo[i]->scale(x);
  } else {
    G4cout << "Histo::scale: WARNING! wrong histogram index " << i << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Histo::print(G4int i)
{
  G4cout<<"### Histogram  "<<i<<"  ###"<<G4endl;
#ifdef G4ANALYSIS_USE
  if(i>=0 && i<nHisto) {
    G4double step = (xmax[i] - xmin[i])/G4double( bins[i]);
    G4double x    =  xmin[i] - step*0.5;
    G4double y, maxX=0, maxY=0;
    G4int    maxJ=0;

    for(G4int j=0; j<bins[i]; j++) {
      x += step;
      y  = histo[i]->binHeight(j);
      if(maxY < y) {maxY = y; maxX = x; maxJ = j;}

      G4cout<<x<<"  "<<y<<G4endl;
    }
    G4cout<<" maxJ  "<<maxJ<<"  maxX  "<<maxX<<"  maxY  "<<maxY<<G4endl;
  }
#endif
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


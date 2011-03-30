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
#include <AIDA/AIDA.h>
#include "HistoMessenger.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::Histo(G4int ver)
{
  verbose    = ver;
  histName   = "histo";
  histType   = "root";
  option     = "";
  nHisto     = 0;
  defaultAct = 1;
  tupleName  = "tree.root";
  tupleId    = "100";
  tupleList  = "";
  ntup       = 0;
  messenger  = 0;

#ifdef G4ANALYSIS_USE
  messenger = new HistoMessenger(this);
  tree = 0;
  af   = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::~Histo()
{
#ifdef G4ANALYSIS_USE
  delete messenger;
  delete af;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::book()
{
#ifdef G4ANALYSIS_USE
  G4cout << "### Histo books " << nHisto << " histograms " << G4endl;
  // Creating the analysis factory
  if(!af) af = AIDA_createAnalysisFactory();
  if(verbose>0) {
    G4cout<<"HIsto books analysis factory ......... "<<G4endl;
  }
  // Creating the tree factory
  AIDA::ITreeFactory* tf = af->createTreeFactory();
  if(verbose>0) {
    G4cout<<"Histo books tree factory ......... "<<G4endl;
  }
  G4String histExt = "";
  char* path = getenv("PHYSLIST");
  if (path) { histExt = "_" + G4String(path); }

  G4String histDir = "";
  path = getenv("HISTODIR");
  if (path) histDir = G4String(path) + "/";

  // Creating a tree mapped to a new hbook file.
  G4String name = histDir + histName + histExt + "." + histType;
  tree = tf->create(name, histType, false, true, option);
  G4cout << "Histo: tree store : " << tree->storeName() << G4endl;
  delete tf;

  // Creating a histogram factory, whose histograms will be handled by the tree
  AIDA::IHistogramFactory* hf = af->createHistogramFactory( *tree );

  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<nHisto; i++) {
    if(active[i]) {
      if(verbose>0) {
	G4cout<<" I am in book: histogram "<< i << " id= " << ids[i] <<G4endl;
      }
      G4String idd;
      if(histType == "root") { idd = "h" +  ids[i]; }
      else                   { idd = ids[i]; }
      histo[i] = hf->createHistogram1D(idd, titles[i], bins[i], xmin[i], xmax[i]);
    } else {
      histo[i] = 0;
    }
  }
  delete hf;
  // Creating a tuple factory, whose tuples will be handled by the tree
  if(tupleList != "") {
    if(verbose>0) {
      G4cout<<"Histo books tuple factory for "<<tupleName <<G4endl;
    }
    AIDA::ITupleFactory* tpf = af->createTupleFactory( *tree );
    ntup = tpf->create(tupleId, tupleName, tupleList);
    delete tpf;
  }
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

void Histo::reset()
{
#ifdef G4ANALYSIS_USE
  delete tree;
  tree = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setFileType(const G4String& nam) 
{
  if(nam == "hbook" || nam == "root" || nam == "aida") histType = nam;
  else if(nam == "XML" || nam == "xml") histType = "aida";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::add1D(const G4String& id, const G4String& name, G4int nb,
                  G4double x1, G4double x2, G4double u)
{
  if(nHisto > 0) {
    for(G4int i=0; i<nHisto; i++) {
      if(ids[i] == id) return;
    }
  }

  if(verbose > 0) 
    G4cout << "New histogram will be booked: #" << id << "  <" << name
           << "  " << nb << "  " << x1 << "  " << x2 << "  " << u
           << G4endl;

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double u)
{
  if(i>=0 && i<nHisto) {
    if(verbose > 0) 
    {
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

void Histo::addTuple(const G4String& w1, const G4String& w2, const G4String& w3)
{
  tupleId = w1;
  tupleName = w2;
  tupleList = w3;

  G4cout<<" addTuple: Id "<<w1<<" Name "<<w2<<" List "<<w3<<G4endl;
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


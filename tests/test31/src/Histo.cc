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

#include <memory> // for the auto_ptr(T>
#include "AIDA/AIDA.h"
#include "HistoMessenger.hh"

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo* Histo::m_instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo* Histo::GetInstance()
{
  if(m_instance == 0){
    static Histo hist; 
    m_instance = &hist;
  }
  return m_instance;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::Histo()
{
  verbose    = 0;
  histName   = "histo.paw";
  histType   = "hbook";
  nHisto     = 0;
  defaultAct = true;
  tupleName  = "tuple.paw";
  tupleId    = "100";
  tupleList  = "";
  ntup       = 0;
  messenger  = 0;

#ifdef G4ANALYSIS_USE
  messenger = new HistoMessenger(this);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::~Histo()
{
#ifdef G4ANALYSIS_USE
  delete messenger;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::clear()
{
#ifdef G4ANALYSIS_USE
  for(G4int i=0; i<nHisto; i++) {
    histo[i]  = 0;
    active[i] = false;
  }
  if(ntup) ntup = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::book()
{
#ifdef G4ANALYSIS_USE
  G4cout << "### Histo books new histograms " << G4endl;
  // Creating the analysis factory
  std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );
  // Creating the tree factory
  std::auto_ptr< AIDA::ITreeFactory > tf( af->createTreeFactory() );

  // Creating a tree mapped to a new file.
  G4String comp = "uncompress";
  if(histType == "xml" || histType == "XML") comp = "compress";
  tree = tf->create(histName,histType,false,true,comp);
  if(tree) 
    G4cout << "Tree store  : " << tree->storeName() << G4endl;
  else
    G4cout << "ERROR: Tree store " << histName  << " is not created!" << G4endl;

  // Creating a histogram factory, whose histograms will be handled by the tree
  std::auto_ptr< AIDA::IHistogramFactory > hf(af->createHistogramFactory( *tree ));

  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<nHisto; i++) {
    if(active[i]) {
      histo[i] = hf->createHistogram1D(ids[i], titles[i], bins[i], xmin[i], xmax[i]);
      if(histo[i] && verbose>0) ListHistogram(i);
    }
  }
  // Creating a tuple factory, whose tuples will be handled by the tree
  if(tupleList != "") {
     std::auto_ptr< AIDA::ITupleFactory > tpf( af->createTupleFactory( *tree ) );
     ntup = tpf->create(tupleId, tupleName, tupleList);
  }
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::save()
{
#ifdef G4ANALYSIS_USE
  // Write histogram file
  tree->commit();
  if(verbose>0) G4cout << "Closing the tree " << G4endl;
  tree->close();
  G4cout << "Histograms and Ntuples are saved" << G4endl;
  clear();
  if(verbose>0) G4cout << "Tree is deleting " << G4endl;
  delete tree;
  if(verbose>0) G4cout << "Tree is deleted " << G4endl;
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Histo::add1D(const G4String& id, const G4String& name, G4int nb, 
                         G4double x1, G4double x2, G4double u)
{
  nHisto++;
  G4double uu = 1.0;
  if(u > 0.0) uu = u;
  if(verbose > 0) {
    G4cout << "New histogram will be booked: id= " << nHisto-1 
           << "  " << id << "  <" << name 
           << ">  " << nb << "  " << x1/uu << "  " << x2/uu << "  " << uu 
           << G4endl;
  }
  active.push_back(defaultAct);
  bins.push_back(nb);
  xmin.push_back(x1);
  xmax.push_back(x2);
  unit.push_back(uu);
  ids.push_back(id);
  titles.push_back(name);
  histo.push_back(0);
  return nHisto-1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double u)
{
  if(i>=0 && i<nHisto) {
    G4double uu = unit[i];
    if(u > 0.0) uu = u;
    if(verbose > 0) {
      G4cout << "Update histogram: #" << i  
             << "  " << nb << "  " << x1/uu << "  " << x2/uu << "  " << uu 
             << G4endl;
    }
    bins[i] = nb;
    xmin[i] = x1;
    xmax[i] = x2;
    unit[i] = uu;
  } else {
    G4cout << "Histo::setHisto1D: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::activate(G4int i, G4bool w)
{
  if(verbose > 0) {
    G4cout << "histogram: id= " << i 
           << "  activation= " << w
           << G4endl;   
  }
  if(i>=0 && i<nHisto) active[i] = w;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::fill(G4int i, G4double x, G4double w)
{
  if(verbose > 1) {
    G4cout << "fill histogram: id= " << i << " at x= " << x 
           << "  weight= " << w << " " << histo[i] 
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
    G4cout << "Scale histogram: id= " << i << " by factor " << x << G4endl;   
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
  if(nam == "hbook" || nam == "root" || nam == "xml" || nam == "XML")
    histType = nam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::ListHistogram(G4int val)
{
  G4int imin = val;
  G4int imax = val;
  if(val < 0 || val >= nHisto) {
    imin = 0;
    imax = nHisto-1;
  }
  for(G4int i=imin; i<=imax; i++) {
    G4cout << "### Histogramm #" << i << " " << ids[i] << " " 
           << titles[i] << G4endl; 
    G4cout << "                       nbins= " << bins[i]
           << " xmin= " << xmin[i]/unit[i]
           << " xmax= " << xmax[i]/unit[i]
           << " unit= " << unit[i]
           << " flag= " << active[i]
           << " " << histo[i]
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Histo::NumberOfBins(G4int id)
{
  G4int n = 0;
  if(id>=0 && id<nHisto) n = bins[id];
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double Histo::MinBin(G4int id)
{
  G4double x = 1.0e-9;
  if(id>=0 && id<nHisto) x = xmin[id];
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4double Histo::MaxBin(G4int id)
{
  G4double x = 1.0e-9;
  if(id>=0 && id<nHisto) x = xmax[id];
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool Histo::IsActive(G4int id)
{
  G4bool yes = false;
  if(id>=0 && id<nHisto) yes = active[id];
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


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

#include "exrdmHisto.hh"
#include "exrdmHistoMessenger.hh"

#ifdef G4ANALYSIS_USE
//#include <memory> // for the auto_ptr(T>
#include <AIDA/AIDA.h>
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmHisto::exrdmHisto()
{
  verbose    = 0;
  //  histName   = "histo.hbook";
  //  histType   = "hbook";
  histName   = "histo.aida";
  histType   = "xml";
  nHisto     = 0;
  nTuple     = 0;
  defaultAct = 1;
  //
  histo.clear();
  ntup.clear();
  active.clear();
  bins.clear();
  xmin.clear();
  xmax.clear();
  unit.clear();
  ids.clear();
  titles.clear();
  tupleName.clear();
  tupleId.clear();
  tupleList.clear();
  messenger  = 0;

  messenger = new exrdmHistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmHisto::~exrdmHisto()
{
#ifdef G4ANALYSIS_USE
  for(G4int i=0; i<nHisto; i++) {
    if(histo[i]) delete histo[i];
  }
  delete messenger;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::book()
{
#ifdef G4ANALYSIS_USE
  G4cout << "### exrdmHisto books " << nHisto << " histograms " << G4endl;
  // Creating the analysis factory
  //  std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );
  AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
  // Creating the tree factory
  //std::auto_ptr< AIDA::ITreeFactory > tf( af->createTreeFactory() );
  AIDA::ITreeFactory* tf = af->createTreeFactory(); 
  // Creating a tree mapped to a new hbook file.

  tree = tf->create(histName,histType,false,true,"uncompress");
  if(tree) 
    G4cout << "Tree store  : " << tree->storeName() << G4endl;
  else
    G4cout << "ERROR: Tree store " << histName  << " is not created!" << G4endl;

  // Creating a histogram factory, whose histograms will be handled by the tree
  //std::auto_ptr< AIDA::IexrdmHistogramFactory > hf(af->createexrdmHistogramFactory( *tree ));
  AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);
  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<nHisto; i++) {
    if(active[i]) {
      histo[i] = hf->createHistogram1D(ids[i], titles[i], bins[i], xmin[i], xmax[i]);
    }
  }
  // Creating a tuple factory, whose tuples will be handled by the tree  
  //std::auto_ptr< AIDA::ITupleFactory > tpf( af->createTupleFactory( *tree ) );
  AIDA::ITupleFactory* tpf =  af->createTupleFactory( *tree );
  for(G4int i=0; i<nTuple; i++) {
    if(tupleList[i] != "") {
      G4cout << "Creating Ntuple: " << tupleName[i] << G4endl;
      ntup[i] = tpf->create(tupleId[i], tupleName[i], tupleList[i]);
    }
  }
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::save()
{
#ifdef G4ANALYSIS_USE
  // Write histogram file
  tree->commit();
  G4cout << "Closing the tree..." << G4endl;
  tree->close();
  G4cout << "exrdmHistograms and Ntuples are saved" << G4endl;
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::add1D(const G4String& id, const G4String& name, G4int nb, 
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::setHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double u)
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
    G4cout << "exrdmHisto::setexrdmHisto1D: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::fillHisto(G4int i, G4double x, G4double w)
{
  if(verbose > 1) {
    G4cout << "fill histogram: #" << i << " at x= " << x 
           << "  weight= " << w
           << G4endl;   
  }
#ifdef G4ANALYSIS_USE  
  if(i>=0 && i<nHisto) {
    histo[i]->fill((float)(x/unit[i]), (float)w);
  } else {
    G4cout << "exrdmHisto::fill: WARNING! wrong histogram index " << i << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::scaleHisto(G4int i, G4double x)
{
  if(verbose > 0) {
    G4cout << "Scale histogram: #" << i << " by factor " << x << G4endl;   
  }
#ifdef G4ANALYSIS_USE  
  if(i>=0 && i<nHisto) {
    histo[i]->scale(x);
  } else {
    G4cout << "exrdmHisto::scale: WARNING! wrong histogram index " << i << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::addTuple(const G4String& w1, const G4String& w2, const G4String& w3)
{
  nTuple++;
  tupleId.push_back(w1);
  tupleName.push_back(w2) ;
  tupleList.push_back(w3) ;
  ntup.push_back(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::fillTuple(G4int i, const G4String& parname, G4double x)
{
  if(verbose > 1) {
    G4cout << "fill tuple # " << i 
	   <<" with  parameter <" << parname << "> = " << x << G4endl; 
  }
#ifdef G4ANALYSIS_USE  
  if(ntup[i]) ntup[i]->fill(ntup[i]->findColumn(parname), (float)x);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::fillTuple(G4int i, const G4String& parname, G4String x)
{
  if(verbose > 1) {
    G4cout << "fill tuple # " << i 
	   <<" with  parameter <" << parname << "> = " << x << G4endl; 
  }
#ifdef G4ANALYSIS_USE  
  if(ntup[i]) ntup[i]->fill(ntup[i]->findColumn(parname), x);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::addRow(G4int i)
{
#ifdef G4ANALYSIS_USE
  if(ntup[i]) ntup[i]->addRow();
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::setFileName(const G4String& nam) 
{
  histName = nam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::setFileType(const G4String& nam) 
{
  histType = nam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


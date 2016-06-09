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
#ifdef G4ANALYSIS_USE
#include <AIDA/AIDA.h>
#endif
//
#ifdef G4ANALYSIS_USE_ROOT
#include "TROOT.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtuple.h"
#endif

#include "exrdmHisto.hh"
#include "exrdmHistoMessenger.hh"
#include "G4ParticleTable.hh"

#include "G4Tokenizer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
exrdmHisto::exrdmHisto()
{
  verbose    = 1;
  histName   = "exrdm";
  //  histType   = "aida";
  histType   = "root";
  nHisto     = 0;
  nTuple     = 0;
  defaultAct = 1;
  //
#ifdef G4ANALYSIS_USE
  aida = 0;
  tree = 0;
#endif

#ifdef G4ANALYSIS_USE_ROOT
  ROOThisto.clear();
  ROOTntup.clear();
  Rarray.clear();
  Rcol.clear();
#endif

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
  tupleListROOT.clear();

  messenger = new exrdmHistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmHisto::~exrdmHisto()
{
#ifdef G4ANALYSIS_USE
  histo.clear();
  ntup.clear();
#endif
#ifdef G4ANALYSIS_USE_ROOT
  //FIXME : G.Barrand : the below is crashy.
  //        In principle the TH are deleted
  //        when doing the TFile::Close !
  //         In fact the hfileROOT should 
  //        be deleted in save(). And I am pretty
  //        sure that the TApplication is not needed.
  //
  // removed by F.Lei
  //  for(G4int i=0; i<nHisto; i++) {
  //   if(ROOThisto[i]) delete ROOThisto[i];
  // }
#endif
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::book()
{
#ifdef G4ANALYSIS_USE
  G4cout << "### exrdmHisto books " << nHisto << " histograms " << G4endl; 
  // Creating the analysis factory
  aida = AIDA_createAnalysisFactory();
  if(!aida) {
    G4cout << "ERROR: can't get AIDA." << G4endl; 
    return;
  }
  // Creating the tree factory
 {AIDA::ITreeFactory* tf = aida->createTreeFactory(); 
  // Creating a tree mapped to a new aida file.
  G4String fileName = histName + "." + histType;
  if (histType == "root") fileName = histName + "_aida." + histType;
  tree = tf->create(fileName,histType,false,true,"");
  delete tf;
  if(!tree) { 
    G4cout << "ERROR: Tree store " << histName  << " is not created!" << G4endl; 
    return;
  }
  G4cout << "Tree store  : " << tree->storeName() << G4endl;}
  // Creating a histogram factory, whose histograms will be handled by the tree
 {AIDA::IHistogramFactory* hf = aida->createHistogramFactory(*tree);
  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<nHisto; i++) {
    if(active[i]) {
      if(verbose>1)
        G4cout<<"book: histogram "<< i << " id= " << ids[i] <<G4endl;
      G4String tit = ids[i];
      if(histType == "root") tit = "h" + ids[i];
      histo[i] = hf->createHistogram1D(tit, titles[i], bins[i], xmin[i], xmax[i]);
    }
  }
  delete hf;
  G4cout << "AIDA histograms are booked" << G4endl;}

  // Creating a tuple factory, whose tuples will be handled by the tree  
 {AIDA::ITupleFactory* tpf =  aida->createTupleFactory( *tree );
  G4cout << "AIDA will book " << nTuple << " ntuples" << G4endl;
  for(G4int i=0; i<nTuple; i++) {
    if(tupleList[i] != "") {
      G4cout << "Creating Ntuple: " << tupleName[i] <<":" <<tupleList[i] <<G4endl;
      ntup[i] = tpf->create(tupleId[i], tupleName[i], tupleList[i],"");
    }
  }
  delete tpf;
  G4cout << "AIDA ntuples are booked" << G4endl;}
#endif

#ifdef G4ANALYSIS_USE_ROOT
//  new TApplication("App", ((int *)0), ((char **)0));
  G4String fileNameROOT = histName + G4String(".root");
  hfileROOT = new TFile(fileNameROOT.c_str() ,"RECREATE","ROOT file for exRDM");
  G4cout << "Root file: " << fileNameROOT << G4endl;
  // Creating an 1-dimensional histograms in the root directory of the tree
  for(G4int i=0; i<nHisto; i++) {
    if(active[i]) {
      G4String id = G4String("h")+ids[i];
      ROOThisto[i] = new TH1D(id, titles[i], bins[i], xmin[i], xmax[i]);
      G4cout << "ROOT Histo " << ids[i] << " " << titles[i] << " booked " << G4endl;
    }
  }
  // Now the ntuples  
  for(G4int i=0; i<nTuple; i++) {
    if(tupleListROOT[i] != "") {
      G4String id = G4String("t")+tupleId[i];
      G4cout << "Creating Ntuple "<<tupleId[i] << " in ROOT file: " 
	     << tupleName[i] << G4endl;
      ROOTntup[i] = new TNtuple(id, tupleName[i], tupleListROOT[i]);
      G4cout << "ROOT Ntuple " << id << " " << tupleName[i] <<" "<< tupleListROOT[i]<< " booked " << G4endl;
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
  G4cout << "Closing the AIDA tree..." << G4endl;
  tree->close();
  G4cout << "Histograms and Ntuples are saved" << G4endl;
  delete tree;
  tree = 0;
  delete aida;
  aida = 0;
  {for(G4int i=0; i<nHisto; i++) histo[i] = 0;}
  {for(G4int i=0; i<nTuple; i++) ntup[i] = 0;}
#endif
#ifdef G4ANALYSIS_USE_ROOT
  G4cout << "ROOT: files writing..." << G4endl;
  hfileROOT->Write();
  G4cout << "ROOT: files closing..." << G4endl;
  hfileROOT->Close();
  //
  // F.Lei added following Guy's suggestion!
  delete hfileROOT;

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
#ifdef G4ANALYSIS_USE
  histo.push_back(0);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  ROOThisto.push_back(0);
#endif
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
    histo[i]->fill(x/unit[i], w);
  } else {
    G4cout << "exrdmHisto::fill: WARNING! wrong AIDA histogram index " << i << G4endl;
  }
#endif
#ifdef G4ANALYSIS_USE_ROOT  
  if(i>=0 && i<nHisto) {
    ROOThisto[i]->Fill(x/unit[i],w);
  } else {
    G4cout << "exrdmHisto::fill: WARNING! wrong ROOT histogram index " << i << G4endl;
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
    G4cout << "exrdmHisto::scale: WARNING! wrong AIDA histogram index " << i << G4endl;
  }
#endif
#ifdef G4ANALYSIS_USE_ROOT  
  if(i>=0 && i<nHisto) {
    ROOThisto[i]->Scale(x);
  } else {
    G4cout << "exrdmHisto::scale: WARNING! wrong ROOT histogram index " << i << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
void exrdmHisto::addTuple(const G4String& w1, const G4String& w2, const G4String& w3 )
#else
void exrdmHisto::addTuple(const G4String& w1, const G4String& w2, const G4String& )
#endif

{
  //G4cout << w1 << " " << w2 << " " << w3 << G4endl;
  nTuple++;
  tupleId.push_back(w1);
  tupleName.push_back(w2) ;
#ifdef G4ANALYSIS_USE
  tupleList.push_back(w3);
  ntup.push_back(0);
#endif

#ifdef G4ANALYSIS_USE_ROOT
  std::vector<float> ar;
  ar.clear();
  for (size_t i = 0; i < 20; i++) ar.push_back(0.);
  Rarray.push_back(ar);
  // convert AIDA header to ROOT header for ntuple
  G4Tokenizer next(w3);
  G4String token = next();
  G4String ROOTList1 = "" ;
  G4int col = 0;
  while ( token != "") {
   token = next();
   if (token == ",") token = next();
   if (token.contains(",")) token.remove(token.size()-1);
   ROOTList1 = ROOTList1 + token + G4String(":");
   col++;
  }
  G4String ROOTList = ROOTList1.substr(0,ROOTList1.length()-2);
//  G4cout << ROOTList << G4endl;
  tupleListROOT.push_back(ROOTList);
  ROOTntup.push_back(0);
  Rcol.push_back(col-1);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::fillTuple(G4int i, const G4String& parname, G4double x)
{
  if(verbose > 1) 
    G4cout << "fill tuple # " << i 
	   <<" with  parameter <" << parname << "> = " << x << G4endl; 
#ifdef G4ANALYSIS_USE
  if(ntup[i]) ntup[i]->fill(ntup[i]->findColumn(parname), x);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::fillTuple(G4int i, G4int col, G4double x)
{
  if(verbose > 1) {
    G4cout << "fill tuple # " << i 
	   <<" in column < " << col << "> = " << x << G4endl; 
  }
#ifdef G4ANALYSIS_USE
  if(ntup[i]) ntup[i]->fill(col,double(x));
#endif

#ifdef G4ANALYSIS_USE_ROOT  
  if(ROOTntup[i]) (Rarray[i])[col] = float(x);
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::fillTuple(G4int i, const G4String& parname, G4String& x)
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
  if(verbose > 1) G4cout << "Added a raw #" << i << " to tuple" << G4endl; 
#ifdef G4ANALYSIS_USE
  if(ntup[i]) ntup[i]->addRow();
#endif

#ifdef G4ANALYSIS_USE_ROOT
  float ar[4];
  for (G4int j=0; j < Rcol[i]; j++) {
//      G4cout << i << " " << Rarray[i][j] << G4endl;
      ar[j] = Rarray[i][j];       
  }  
  if(ROOTntup[i]) ROOTntup[i]->Fill(ar);
#endif

} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::setFileName(const G4String& nam) 
{
  histName = nam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4String& exrdmHisto::getFileName() const
{
  return histName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmHisto::setFileType(const G4String& nam) 
{
  histType = nam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4String& exrdmHisto::FileType() const
{
  return histType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


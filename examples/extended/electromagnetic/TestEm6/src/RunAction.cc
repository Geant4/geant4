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
//
// $Id: RunAction.cc,v 1.6 2004-03-15 11:23:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"

#ifdef USE_AIDA
 #include "AIDA/AIDA.h"
#endif

#ifdef USE_ROOT
 #include "TFile.h"
 #include "TH1F.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
 cleanHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::bookHisto()
{
#ifdef USE_AIDA
 // Creating the analysis factory
 AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
 AIDA::ITreeFactory* tf = af->createTreeFactory();
 
 // Creating a tree mapped to an hbook file.
 G4bool readOnly  = false;
 G4bool createNew = true;
 tree = tf->create("testem6.paw", "hbook", readOnly, createNew);

 // Creating a histogram factory, whose histograms will be handled by the tree
 AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);

 // Creating histograms
 histo[0] = hf->createHistogram1D("1","1/(1+(theta+[g]+)**2)",100, 0 ,1.);
 histo[1] = hf->createHistogram1D("2","log10(theta+ [g]+)",   100,-3.,1.);
 histo[2] = hf->createHistogram1D("3","log10(theta- [g]-)",   100,-3.,1.);
 histo[3] = hf->createHistogram1D("4","log10(theta+ [g]+ -theta- [g]-)",100,-3.,1.);
 histo[4] = hf->createHistogram1D("5","xPlus" ,100,0.,1.);
 histo[5] = hf->createHistogram1D("6","xMinus",100,0.,1.);
  
 delete hf;
 delete tf;
 delete af;     
#endif

#ifdef USE_ROOT
 // Create a ROOT file
 tree = new TFile("testem7.root","recreate");
 
 // Create the histogram
 histo[0] = new TH1F("1","1/(1+(theta+[g]+)**2)",100, 0 ,1.);
 histo[1] = new TH1F("2","log10(theta+ [g]+)",   100,-3.,1.);
 histo[2] = new TH1F("3","log10(theta- [g]-)",   100,-3.,1.);
 histo[3] = new TH1F("4","log10(theta+ [g]+ -theta- [g]-)",100,-3.,1.);
 histo[4] = new TH1F("5","xPlus" ,100,0.,1.);
 histo[5] = new TH1F("6","xMinus",100,0.,1.); 
#endif
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::cleanHisto()
{
#ifdef USE_AIDA
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)  
  delete tree;
#endif
  
#ifdef USE_ROOT
  tree->Write();        // Writing the histograms to the file
  tree->Close();        // and closing the file  
  delete tree;
#endif     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();
     
  //histograms
  //
  if (aRun->GetRunID() == 0) bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

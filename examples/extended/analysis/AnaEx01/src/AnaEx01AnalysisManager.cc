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
//
// $Id: AnaEx01AnalysisManager.cc,v 1.16.4.1 2009/08/11 09:50:10 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-02 $
//
// 
// Guy Barrand 14th Septembre 2000.

#ifdef G4ANALYSIS_USE

#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/ITreeFactory.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITuple.h>

#include "AnaEx01CalorHit.hh"
#include "AnaEx01AnalysisManager.hh"

AnaEx01AnalysisManager::AnaEx01AnalysisManager(AIDA::IAnalysisFactory* aAIDA)
:fCalorimeterCollID(-1)
,fAIDA(aAIDA)
,fTree(0)
,fEAbs(0)
,fLAbs(0)
,fEGap(0)
,fLGap(0)
,fTuple(0)
{
  // Could fail if no AIDA implementation found :
  if(!fAIDA) {
    G4cout << "AIDA analysis factory not found." << G4endl;
    return;
  }

  AIDA::ITreeFactory* treeFactory = fAIDA->createTreeFactory();
  if(!treeFactory) return;

  // Create a tree-like container to handle histograms.
  // This tree is associated to a AnaEx01.<format> file.  

  std::string h_name_EAbs("EAbs");
  std::string h_name_LAbs("LAbs");
  std::string h_name_EGap("EGap");
  std::string h_name_LGap("LGap");
  std::string t_name("AnaEx01");

  // File format :
  std::string format("xml");
  //std::string format("hbook");
  //std::string format("root");

  std::string file("AnaEx01");
  std::string ext;
  ext = "."+format;

  std::string opts;
  if(format=="hbook") {
    opts = "compress=no";

    h_name_EAbs = "1";
    h_name_LAbs = "2";
    h_name_EGap = "3";
    h_name_LGap = "4";
    t_name = "101";

  } else if(format=="root") {
    opts = "compress=yes";

  } else if(format=="xml") {
    ext = ".aida";
    opts = "compress=no";

  } else {
    G4cout << "storage format \"" << format << "\""
           << " not handled in this example."
           << G4endl;
    return;
  }

  file += ext;
  fTree = treeFactory->create(file,format,false,true,opts);

  // Factories are not "managed" by an AIDA analysis system.
  // They must be deleted by the AIDA user code.
  delete treeFactory; 

  if(!fTree) {
    G4cout << "can't create tree associated to file \"" << file << "\"."
           << G4endl;
    return;
  }

  fTree->mkdir("histograms");
  fTree->cd("histograms");
      
  // Create an histo factory that will create histo in the tree :
  AIDA::IHistogramFactory* histoFactory = 
    fAIDA->createHistogramFactory(*fTree);
  if(histoFactory) {
    fEAbs = histoFactory->createHistogram1D(h_name_EAbs,"EAbs",100,0,100);
    if(!fEAbs) G4cout << "can't create histo EAbs." << G4endl;
    fLAbs = histoFactory->createHistogram1D(h_name_LAbs,"LAbs",100,0,100);
    if(!fLAbs) G4cout << "can't create histo LAbs." << G4endl;
    fEGap = histoFactory->createHistogram1D(h_name_EGap,"EGap",100,0,10);
    if(!fEGap) G4cout << "can't create histo EGap." << G4endl;
    fLGap = histoFactory->createHistogram1D(h_name_LGap,"LGap",100,0,100);
    if(!fLGap) G4cout << "can't create histo LGap." << G4endl;
    delete histoFactory;
  }
    
  fTree->cd("..");
  fTree->mkdir("tuples");
  fTree->cd("tuples");
    
  // Get a tuple factory :
  AIDA::ITupleFactory* tupleFactory = fAIDA->createTupleFactory(*fTree);
  if(tupleFactory) {
    
    // Create a tuple :
    fTuple = tupleFactory->create(t_name,"AnaEx01",
      "double EAbs,double LAbs,double EGap,double LGap");
    if(!fTuple) G4cout << "can't create tuple." << G4endl;
    
    delete tupleFactory;
  }

  fTree->cd("..");
}
AnaEx01AnalysisManager::~AnaEx01AnalysisManager() {
}

void AnaEx01AnalysisManager::BeginOfRun(const G4Run* aRun){
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

void AnaEx01AnalysisManager::EndOfRun(const G4Run*){
  if(fTree) fTree->commit();
  if(fEAbs) {
    G4cout << "Histo : EAbs : mean " << fEAbs->mean() << " rms : " << fEAbs->rms() << G4endl;
    G4cout << "Histo : LAbs : mean " << fLAbs->mean() << " rms : " << fLAbs->rms() << G4endl;
    G4cout << "Histo : EGap : mean " << fEGap->mean() << " rms : " << fEGap->rms() << G4endl;
    G4cout << "Histo : LGap : mean " << fLGap->mean() << " rms : " << fLGap->rms() << G4endl;
  }
}

void AnaEx01AnalysisManager::BeginOfEvent(const G4Event*){
  if(fCalorimeterCollID==-1) {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    fCalorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 
}

void AnaEx01AnalysisManager::EndOfEvent(const G4Event* aEvent){
  if(!fEAbs) return; // No histo booked !
  if(!fTuple) return; // No tuple booked !

  //G4int evtNb = aEvent->GetEventID();

  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  AnaEx01CalorHitsCollection* CHC = 
    HCE ? (AnaEx01CalorHitsCollection*)(HCE->GetHC(fCalorimeterCollID)) : 0;

  if(CHC) {
    G4int n_hit = CHC->entries();
    for (G4int i=0;i<n_hit;i++) {
      G4double EAbs = (*CHC)[i]->GetEdepAbs();
      G4double LAbs = (*CHC)[i]->GetTrakAbs();
      G4double EGap = (*CHC)[i]->GetEdepGap();
      G4double LGap = (*CHC)[i]->GetTrakGap();
      fEAbs->fill(EAbs);
      fLAbs->fill(LAbs);
      fEGap->fill(EGap);
      fLGap->fill(LGap);

      fTuple->fill(0,EAbs);
      fTuple->fill(1,LAbs);
      fTuple->fill(2,EGap);
      fTuple->fill(3,LGap);
      fTuple->addRow();
    }
  }	
  
}
void AnaEx01AnalysisManager::Step(const G4Step*){}

#endif

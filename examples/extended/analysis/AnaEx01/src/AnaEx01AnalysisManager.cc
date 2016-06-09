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
// $Id: AnaEx01AnalysisManager.cc,v 1.15 2005/10/25 18:19:56 gbarrand Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
  // This tree is associated to a AnaEx01.aida file.  
  //std::string opts = "compress=yes";
  std::string opts = "compress=no";
  fTree = treeFactory->create("AnaEx01.aida","xml",false,true,opts);
  //std::string opts = "export=root";
  //fTree = treeFactory->create("AnaEx01.root","ROOT",false,true,opts);

  // Factories are not "managed" by an AIDA analysis system.
  // They must be deleted by the AIDA user code.
  delete treeFactory; 

  if(!fTree) return;

  fTree->mkdir("histograms");
  fTree->cd("histograms");
      
  // Create an histo factory that will create histo in the tree :
  AIDA::IHistogramFactory* histoFactory = 
    fAIDA->createHistogramFactory(*fTree);
  if(histoFactory) {
    fEAbs = histoFactory->createHistogram1D("EAbs",100,0,100);
    fLAbs = histoFactory->createHistogram1D("LAbs",100,0,100);
    fEGap = histoFactory->createHistogram1D("EGap",100,0,10);
    fLGap = histoFactory->createHistogram1D("LGap",100,0,100);
    delete histoFactory;
  }
    
  fTree->cd("..");
  fTree->mkdir("tuples");
  fTree->cd("tuples");
    
  // Get a tuple factory :
  AIDA::ITupleFactory* tupleFactory = fAIDA->createTupleFactory(*fTree);
  if(tupleFactory) {
    
    // Create a tuple :
    fTuple = tupleFactory->create("AnaEx01","AnaEx01",
      "double EAbs,double LAbs,double EGap,double LGap");
    
    delete tupleFactory;
  }

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

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
// CaTS (Calorimetry and Tracking Simulation)
//
// Authors: Hans Wenzel and Soon Yung Jun
//          (Fermi National Accelerator Laboratory)
//
// History:
// October  18th, 2021 : first implementation
// November 10th, 2021 : implement writing one file per worker thread
//                       which are then merged
// ********************************************************************
//
/// \file RootIO.cc
/// \brief Implementation of the CaTS::RootIO class

#ifdef WITH_ROOT
// Project headers
#include "RootIO.hh"
#include "ConfigurationManager.hh"
#include "Event.hh"
#include "PhotonSD.hh"
#include "PhotonHit.hh"
#include "InteractionHit.hh"
#include "lArTPCHit.hh"
#include "TrackerHit.hh"
#include "MscHit.hh"
#include "CalorimeterHit.hh"
#include "DRCalorimeterHit.hh"
// Geant4 headers
#include "G4Event.hh"
#include "G4ios.hh"
#include "G4Threading.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHit.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
// Root headers
#include "TBranch.h"
#include "TFile.h"
#include "TObject.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

// G4Mutex in a file scope
namespace
{
  // Mutex to lock RootIO constructor
  G4Mutex RootIOMutex = G4MUTEX_INITIALIZER;
} // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RootIO::RootIO()
{
  G4AutoLock lock(&RootIOMutex);

  ROOT::EnableThreadSafety();
  gSystem->Load("libCaTSClassesDict");

  fFileName = ConfigurationManager::getInstance()->getfname() + ".root";

  if (G4Threading::IsWorkerThread())
  {
    fFileName += std::to_string(G4Threading::G4GetThreadId());
  }

  if (G4Threading::IsWorkerThread() ||
      !G4Threading::IsMultithreadedApplication())
  {
    G4cout << "RootIO:: Opening File: " << fFileName << G4endl;
    fFile = TFile::Open(fFileName.c_str(), "RECREATE");
    fTree = new TTree("Events", "ROOT tree containing Hit collections");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RootIO* RootIO::GetInstance()
{
  static G4ThreadLocal RootIO instance;
  return &instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::Write(const G4Event* event)
{
  G4bool verbose = ConfigurationManager::getInstance()->isEnable_verbose();

  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if (HCE == nullptr)
  {
    return;
  }

  // Record hits from sensitive detectors into the CaTS event data
  Event* CaTSEvt = new Event();
  CaTSEvt->SetEventNr(event->GetEventID());

  if (verbose)
  {
    G4cout << "Number of collections:  " << HCE->GetNumberOfCollections()
           << G4endl;
  }

  for (int i = 0; i < HCE->GetNumberOfCollections(); i++)
  {
    G4VHitsCollection* hc = HCE->GetHC(i);
    std::vector<G4String> y = Split(hc->GetName(), '_');
    G4String Classname = y[1];

    if (verbose)
    {
      G4cout << "Classname: " << Classname << " Collection size: "
	     <<  hc->GetSize() << G4endl;
    }

    if (Classname == "lArTPC")
    {
      AddHits<lArTPCHit>(hc, CaTSEvt);
    }
    else if (Classname == "PhotonDetector")
    {
      AddHits<PhotonHit>(hc, CaTSEvt);
    }
    else if (Classname == "Target")
    {
      AddHits<InteractionHit>(hc, CaTSEvt);
    }
    else if (Classname == "Tracker")
    {
      AddHits<TrackerHit>(hc, CaTSEvt);
    }
    else if (Classname == "Msc")
    {
      AddHits<MscHit>(hc, CaTSEvt);
    }
    else if (Classname == "Calorimeter")
    {
      AddHits<CalorimeterHit>(hc, CaTSEvt);
    }
    else if (Classname == "DRCalorimeter")
    {
      AddHits<DRCalorimeterHit>(hc, CaTSEvt);
    }
    else
    {
      G4cout << "SD type: " << Classname << " unknown" << G4endl;
    }
  }

  if (verbose)
  {
    G4cout << "writing Event: " << CaTSEvt->GetEventNumber() << G4endl;
  }

  G4AutoLock lock(&RootIOMutex);
  if (!fEvtBranch)
  {
    Int_t bufsize = 64000;
    fEvtBranch = fTree->Branch("event.", &CaTSEvt, bufsize, 0);
  }
  else
  {
    fEvtBranch->SetAddress(&CaTSEvt);
  }

  fTree->Fill();
  fFile->Write("", TObject::kOverwrite);

  CaTSEvt->Reset();
  delete CaTSEvt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::Close()
{
  G4cout << " Closing File: " << fFileName << G4endl;
  fFile->Close();
  fEvtBranch = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::Merge()
{
  if (G4Threading::IsMasterThread())
  {
    auto nthreads = G4RunManager::GetRunManager()->GetNumberOfThreads();
    G4cout << "RootIO::Merging: " << nthreads << " threads" << G4endl;

    std::vector<TFile*> files;
    std::vector<TTree*> trees;
    TList* list = new TList;

    for (int i = 0; i < nthreads; i++)
    {
      G4String fileName = fFileName + std::to_string(i);
      files.push_back(TFile::Open(fileName.c_str()));
      trees.push_back((TTree*) (files[i]->Get("Events")));
      list->Add(trees[i]);

      if (i == nthreads - 1)
      {
        auto* file = TFile::Open(fFileName.c_str(), "RECREATE");

        auto tree = TTree::MergeTrees(list);
        tree->SetName("Events");

        Event* event = new Event();
        tree->SetBranchAddress("event.", &event);
        TBranch* branch = tree->GetBranch("event.");
        int nevents = branch->GetEntries();
        G4cout << "Number of Events:  " << nevents << G4endl;

        tree->Write();
        file->Close();
      }
      // Delete the merged file
      if(std::remove(fileName.c_str()) !=0)
      {
	G4cout << "Error deleting file: " << fileName << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4String> RootIO::Split(const G4String& s, char delim)
{
    std::stringstream ss(s);
    G4String item;
    std::vector<G4String> elems;

    while (std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

#endif

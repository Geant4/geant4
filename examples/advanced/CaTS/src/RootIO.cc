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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October  18th, 2021 : first implementation
//   November 10th, 2021 : implement writing one file
//                         per worker thread which are then merged
//
// ********************************************************************
//
/// \file RootIO.cc
/// \brief Implementation of the CaTS::RootIO class

#ifdef WITH_ROOT
// Project headers
#include "RootIO.hh"
#include "ConfigurationManager.hh"
#include "Event.hh"
// Geant4 headers
#include <G4String.hh>
#include <G4ios.hh>
#include <G4Threading.hh>
#include <G4RunManager.hh>
// Root headers
#include "TBranch.h"
#include "TFile.h"
#include "TObject.h"
#include "TSystem.h"
#include "TTree.h"
#include "TROOT.h"
// c++ headers
#include <string>
#include <stdio.h>

G4ThreadLocal RootIO* RootIO::fgInstance = nullptr;
// mutex in a file scope
namespace {
    // Mutex to lock RootIO constructor
    G4Mutex RootIOMutex = G4MUTEX_INITIALIZER;
} // namespace

RootIO::RootIO() {
    G4AutoLock lock(&RootIOMutex);
    TSystem ts;
    gSystem->Load("libCaTSClassesDict");
    if (G4Threading::IsMultithreadedApplication()) {
        if (G4Threading::IsWorkerThread()) {
            G4String fFilename = ConfigurationManager::getInstance()->getfname() +
                    "_Thread_" +
                    std::to_string(G4Threading::G4GetThreadId()) + ".root";
            G4cout << "RootIO:: Opening File: " << fFilename << G4endl;
            fFile = new TFile(fFilename.c_str(), "RECREATE");
            TTree::SetMaxTreeSize(1000 * Long64_t(2000000000));
            // Create a ROOT Tree and one superbranch
            ftree = new TTree("Events", "ROOT tree containing Hit collections");
            ftree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
            Int_t branchStyle = 1;
            TTree::SetBranchStyle(branchStyle);
        }
    } else {
        G4String fFilename = ConfigurationManager::getInstance()->getfname() + ".root";
        G4cout << "RootIO:: Opening File: " << fFilename << G4endl;
        fFile = new TFile(fFilename.c_str(), "RECREATE");
        TTree::SetMaxTreeSize(1000 * Long64_t(2000000000));
        // Create a ROOT Tree and one superbranch
        ftree = new TTree("Events", "ROOT tree containing Hit collections");
        ftree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
        Int_t branchStyle = 1;
        TTree::SetBranchStyle(branchStyle);
    }
}

RootIO* RootIO::GetInstance() {
    if (fgInstance == nullptr) {
        static G4ThreadLocalSingleton<RootIO> inst;
        fgInstance = inst.Instance();
    }
    return fgInstance;
}

void RootIO::Write(Event* fevent) {
    if (ConfigurationManager::getInstance()->isEnable_verbose())
        G4cout << "writing Event: " << fevent->GetEventNumber() << G4endl;
    if ((fevent->GetEventNumber()) % 1000 == 0)
        G4cout << "writing Event: " << fevent->GetEventNumber() << G4endl;
    if (!fevtinitialized) {
        Int_t bufsize = 64000;
        fevtbranch = ftree->Branch("event.", &fevent, bufsize, 0);
        fevtbranch->SetAutoDelete(kFALSE);
        fevtinitialized = true;
    }
    fFile = ftree->GetCurrentFile(); // just in case we switched to a new file
    fnb += ftree->Fill();
    fFile->Write("", TObject::kOverwrite);
}

void RootIO::Close() {
    G4cout << " Closing File: "  << G4endl;
    fFile = ftree->GetCurrentFile();
    fFile->Close();
}

void RootIO::Merge() {
    if (G4Threading::IsMasterThread()) {
        unsigned int nthreads = G4RunManager::GetRunManager()->GetNumberOfThreads();
        G4cout << "RootIO::Merging: " << nthreads << " threads" << G4endl;
        G4String filename = ConfigurationManager::getInstance()->getfname();
        G4String tmpfn;
        std::vector<TFile*> filevec;
        std::vector<Event*> evtvec;
        std::vector<TTree*> treevec;
        std::vector<TBranch*> branchvec;
        TList* list = new TList;
        TTree* newtree;
        for (unsigned int i = 0; i < nthreads; i++) {
            tmpfn = filename + "_Thread_" + std::to_string(i) + ".root";
            filevec.push_back(new TFile(tmpfn.c_str()));
            treevec.push_back((TTree*) (filevec[i]->Get("Events")));
            list->Add(treevec[i]);
            if (i == nthreads - 1) {
                G4String tmpfn2 = filename + "_Merged.root";
                TFile* fm = new TFile(tmpfn2.c_str(), "RECREATE");
                newtree = TTree::MergeTrees(list);
                newtree->SetName("Events");
                Event* eventm = new Event();
                newtree->SetBranchAddress("event.", &eventm);
                TBranch* fevtbranchm = newtree->GetBranch("event.");
                int neventm = fevtbranchm->GetEntries();
                G4cout << "Nr. of Events:  " << neventm << G4endl;
                newtree->Write();
                fm->Close();
            }
        }
        for (unsigned int i = 0; i < nthreads; i++) {
            tmpfn = filename + "_Thread_" + std::to_string(i) + ".root";
            if (remove(tmpfn.c_str()) != 0) {
                G4cout << "Error deleting file: " << tmpfn << G4endl;
            } else {
                G4cout << "File: " << tmpfn << " successfully deleted" << G4endl;
            }
        }
    }
}
#endif

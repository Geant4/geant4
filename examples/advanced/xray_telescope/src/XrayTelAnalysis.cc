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
//
// Author:  A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copied from his UserAnalyser class)
//
// History:
// -----------
// 19 Mar 2013   LP         Migrated to G4tools
//  8 Nov 2002   GS         Migration to AIDA 3
//  7 Nov 2001   MGP        Implementation
//
// -------------------------------------------------------------------

#include <fstream>
#include <iomanip>
#include "G4AutoLock.hh"
#include "XrayTelAnalysis.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"

XrayTelAnalysis* XrayTelAnalysis::instance = 0;

namespace { 
  //Mutex to acquire access to singleton instance check/creation
  G4Mutex instanceMutex = G4MUTEX_INITIALIZER;
  //Mutex to acquire accss to histograms creation/access
  //It is also used to control all operations related to histos 
  //File writing and check analysis
  G4Mutex dataManipulationMutex = G4MUTEX_INITIALIZER;
}

XrayTelAnalysis::XrayTelAnalysis() : 
  asciiFile(0),nEnteringTracks(0), totEnteringEnergy(0)
{
  histFileName = "xraytel";


  asciiFileName="xraytel.out";
  asciiFile = new std::ofstream(asciiFileName);

  if(asciiFile->is_open()) 
    (*asciiFile) << "Energy (keV)  x (mm)    y (mm)    z (mm)" << G4endl << G4endl;  
}

XrayTelAnalysis::~XrayTelAnalysis()
{
  if (asciiFile)
    delete asciiFile;
  if (nEnteringTracks)
    delete nEnteringTracks;
  if (totEnteringEnergy)
    delete totEnteringEnergy;
}


XrayTelAnalysis* XrayTelAnalysis::getInstance()
{
  G4AutoLock l(&instanceMutex);
  if (instance == 0) 
    instance = new XrayTelAnalysis;
  return instance;
}


void XrayTelAnalysis::book(G4bool isMaster)
{
  G4AutoLock l(&dataManipulationMutex);

  //reset counters: do be done only once, by the master
  if (isMaster)
    {
      if (nEnteringTracks)    
	{
	  delete nEnteringTracks;
	  nEnteringTracks = 0;
	}
      nEnteringTracks = new std::map<G4int,G4int>;
  
      if (totEnteringEnergy)
	{
	  delete totEnteringEnergy;
	  totEnteringEnergy = 0;
	}
      totEnteringEnergy = new std::map<G4int,G4double>;
    }

  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Open an output file: it is done in master and threads. The 
  // printout is done only by the master, for tidyness
  if (isMaster)
    G4cout << "Opening output file " << histFileName << " ... ";
  man->OpenFile(histFileName);
  man->SetFirstHistoId(1);
  if (isMaster)
    G4cout << " done" << G4endl;

  // Book 1D histograms
  man->CreateH1("h1","Energy, all /keV",  100,0.,100.);
  man->CreateH1("h2","Energy, entering detector /keV", 500,0.,500.);

  // Book 2D histograms (notice: the numbering is independent)
  man->CreateH2("d1","y-z, all /mm", 100,-500.,500.,100,-500.,500.); 
  man->CreateH2("d2","y-z, entering detector /mm", 200,-50.,50.,200,-50.,50.);
  
  // Book ntuples
  man->CreateNtuple("tree", "Track ntuple");
  man->CreateNtupleDColumn("energy");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleDColumn("dirx");
  man->CreateNtupleDColumn("diry");
  man->CreateNtupleDColumn("dirz");
  man->FinishNtuple();
}


void XrayTelAnalysis::finish(G4bool isMaster)
{
  G4AutoLock l(&dataManipulationMutex);
  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  // Complete clean-up
  delete G4AnalysisManager::Instance();

  if (!isMaster)
    return;

  //only master performs these operations
  asciiFile->close();
 
  //Sequential run: output results!
  if (nEnteringTracks->count(-2))
    {
      G4cout << "End of Run summary (sequential)" << G4endl << G4endl;
      G4cout << "Total Entering Detector : " << nEnteringTracks->find(-2)->second  
	     << G4endl;
      G4cout << "Total Entering Detector Energy : " 
	     << (totEnteringEnergy->find(-2)->second)/MeV  
	     << " MeV"
	     << G4endl;
      return;
    }
  

  //MT run: sum results
 

  //MT build, but sequential run  
  if (nEnteringTracks->count(-1))
    {
      G4cout << "End of Run summary (sequential with MT build)" << G4endl << G4endl;
      G4cout << "Total Entering Detector : " << nEnteringTracks->find(-1)->second  
	     << G4endl;
      G4cout << "Total Entering Detector Energy : " 
	     << (totEnteringEnergy->find(-1)->second)/MeV  
	     << " MeV"
	     << G4endl;
      G4cout << "##########################################" << G4endl;
      return;
    }

  G4bool loopAgain = true;
  G4int totEntries = 0;
  G4double totEnergy = 0.;

  G4cout << "##########################################" << G4endl;
  for (size_t i=0; loopAgain; i++)
    {
      //ok, this thread was found
      if (nEnteringTracks->count(i))
	{
	  G4cout << "End of Run summary (thread= " << i << ")" << G4endl;
	  G4int part = nEnteringTracks->find(i)->second;
	  G4cout << "Total Entering Detector : " << part << G4endl;
	  G4double ene = totEnteringEnergy->find(i)->second;
	  G4cout << "Total Entering Detector Energy : " 
		 << ene/MeV  
		 << " MeV" << G4endl;
	  totEntries += part;
	  totEnergy += ene;
	  G4cout << "##########################################" << G4endl;
	}
      else
	loopAgain = false;
    }
  //Report global numbers
  if (totEntries)
    {
      G4cout << "End of Run summary (MT) global" << G4endl << G4endl;
      G4cout << "Total Entering Detector : " << totEntries << G4endl;
      G4cout << "Total Entering Detector Energy : " 
	     << totEnergy/MeV  
	     << " MeV" << G4endl;
      G4cout << "##########################################" << G4endl;
    }
}

void XrayTelAnalysis::analyseStepping(const G4Track& track, G4bool entering)
{
  G4AutoLock l(&dataManipulationMutex);
  eKin = track.GetKineticEnergy()/keV;
  G4ThreeVector pos = track.GetPosition()/mm;
  y = pos.y();
  z = pos.z();
  G4ThreeVector dir = track.GetMomentumDirection();
  dirX = dir.x();
  dirY = dir.y();
  dirZ = dir.z();

  // Fill histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(1,eKin);
  man->FillH2(1,y,z);
  
  // Fill histograms and ntuple, tracks entering the detector
  if (entering) {
    // Fill and plot histograms
    man->FillH1(2,eKin);
    man->FillH2(2,y,z);

    man->FillNtupleDColumn(0,eKin);
    man->FillNtupleDColumn(1,x);
    man->FillNtupleDColumn(2,y);
    man->FillNtupleDColumn(3,z);
    man->FillNtupleDColumn(4,dirX);
    man->FillNtupleDColumn(5,dirY);
    man->FillNtupleDColumn(6,dirZ);
    man->AddNtupleRow();  
  }


  // Write to file
  if (entering) {
    if(asciiFile->is_open()) {
      (*asciiFile) << std::setiosflags(std::ios::fixed)
		   << std::setprecision(3)
		   << std::setiosflags(std::ios::right)
		   << std::setw(10);
      (*asciiFile) << eKin;
      (*asciiFile) << std::setiosflags(std::ios::fixed)
		   << std::setprecision(3)
		   << std::setiosflags(std::ios::right)
		   << std::setw(10);
      (*asciiFile) << x;
      (*asciiFile) << std::setiosflags(std::ios::fixed)
		   << std::setprecision(3)
		   << std::setiosflags(std::ios::right)
		   << std::setw(10);
      (*asciiFile) << y;
      (*asciiFile) << std::setiosflags(std::ios::fixed)
		   << std::setprecision(3)
		   << std::setiosflags(std::ios::right)
		   << std::setw(10);
      (*asciiFile) << z
		   << G4endl;
    }
  }
}

void XrayTelAnalysis::Update(G4double energy,G4int threadID)
{
  G4AutoLock l(&dataManipulationMutex);
  //It already exists: increase the counter
  if (nEnteringTracks->count(threadID))
    {
      (nEnteringTracks->find(threadID)->second)++;     
    }
  else //enter a new one
    {
      G4int tracks = 1;
      nEnteringTracks->insert(std::make_pair(threadID,tracks));
    }

  //It already exists: increase the counter
  if (totEnteringEnergy->count(threadID))
    (totEnteringEnergy->find(threadID)->second) += energy;
  else //enter a new one    
    totEnteringEnergy->insert(std::make_pair(threadID,energy));
    
  return;
}


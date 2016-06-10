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
// $Id: XrayTelAnalysis.cc 68710 2013-04-05 09:04:21Z gcosmo $
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

#include "XrayTelAnalysis.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"

XrayTelAnalysis* XrayTelAnalysis::instance = 0;

XrayTelAnalysis::XrayTelAnalysis() : 
  asciiFile(0)
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
}


XrayTelAnalysis* XrayTelAnalysis::getInstance()
{
  if (instance == 0) instance = new XrayTelAnalysis;
  return instance;
}


void XrayTelAnalysis::book()
{
  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Open an output file
  G4cout << "Opening output file " << histFileName << " ... ";
  man->OpenFile(histFileName);
  man->SetFirstHistoId(1);
  G4cout << " done" << G4endl;

  // Book 1D histograms
  man->CreateH1("1","Energy, all /keV",  100,0.,100.);
  man->CreateH1("2","Energy, entering detector /keV", 500,0.,500.);

  // Book 2D histograms (notice: the numbering is independent)
  man->CreateH2("1","y-z, all /mm", 100,-500.,500.,100,-500.,500.); 
  man->CreateH2("2","y-z, entering detector /mm", 200,-50.,50.,200,-50.,50.);
  
  // Book ntuples
  man->CreateNtuple("10", "Track ntuple");
  man->CreateNtupleDColumn("energy");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleDColumn("dirx");
  man->CreateNtupleDColumn("diry");
  man->CreateNtupleDColumn("dirz");
  man->FinishNtuple();
}


void XrayTelAnalysis::finish()
{
  asciiFile->close();

  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  // Complete clean-up
  delete G4AnalysisManager::Instance();
}

void XrayTelAnalysis::analyseStepping(const G4Track& track, G4bool entering)
{
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


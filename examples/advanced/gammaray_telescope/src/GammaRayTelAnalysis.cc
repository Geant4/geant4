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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayAnalysisManager  ------
//           by R.Giannitrapani, F.Longo & G.Santin (03 dic 2000)
//
// 03.04.2013 F.Longo/L.Pandola
// - migrated to G4tools
//
// 29.05.2003 F.Longo 
// - anaphe 5.0.5 compliant
//
// 18.06.2002 R.Giannitrapani, F.Longo & G.Santin
// - new release for Anaphe 4.0.3
//
// 07.12.2001 A.Pfeiffer
// - integrated Guy's addition of the ntuple
//
// 06.12.2001 A.Pfeiffer
// - updating to new design (singleton)
//
// 22.11.2001 G.Barrand
// - Adaptation to AIDA
//
// ************************************************************
#include <fstream>
#include <iomanip>

#include "G4RunManager.hh" 

#include "GammaRayTelAnalysis.hh"
#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelAnalysisMessenger.hh"

GammaRayTelAnalysis* GammaRayTelAnalysis::instance = 0;

//-------------------------------------------------------------------------------- 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysis::GammaRayTelAnalysis()
  :GammaRayTelDetector(0),histo2DMode("strip")
{
  GammaRayTelDetector =
    static_cast<const GammaRayTelDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // Define the messenger and the analysis system
  analysisMessenger = new GammaRayTelAnalysisMessenger(this);  
  histoFileName = "gammaraytel";
}


GammaRayTelAnalysis::~GammaRayTelAnalysis() {
  Finish();
  // Complete clean-up
  delete G4AnalysisManager::Instance();
}


void GammaRayTelAnalysis::Init()
{;}                       

void GammaRayTelAnalysis::Finish()
{
  delete analysisMessenger;
  analysisMessenger = 0;
}             

GammaRayTelAnalysis* GammaRayTelAnalysis::getInstance()
{
  if (instance == 0) instance = new GammaRayTelAnalysis();
  return instance;
}

// This function fill the 2d histogram of the XZ positions
void GammaRayTelAnalysis::InsertPositionXZ(double x, double z)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH2(1,x,z);
}

// This function fill the 2d histogram of the YZ positions
void GammaRayTelAnalysis::InsertPositionYZ(double y, double z)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH2(2,y,z);
}

// This function fill the 1d histogram of the energy released in the last Si plane
void GammaRayTelAnalysis::InsertEnergy(double en)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(1,en);
}

// This function fill the 1d histogram of the hits distribution along the TKR planes
void GammaRayTelAnalysis::InsertHits(int nplane)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(2,nplane);
}

void GammaRayTelAnalysis::setNtuple(float E, float p, float x, 
				    float y, float z)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleDColumn(0,E);
  man->FillNtupleDColumn(1,p);
  man->FillNtupleDColumn(2,x);
  man->FillNtupleDColumn(3,y);
  man->FillNtupleDColumn(4,z);
  man->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/

void GammaRayTelAnalysis::BeginOfRun() 
{ 
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Open an output file

  G4cout << "Opening output file " << histoFileName << " ... ";
  man->OpenFile(histoFileName);
  man->SetFirstHistoId(1);
  G4cout << " done" << G4endl;
  

  int Nplane = GammaRayTelDetector->GetNbOfTKRLayers();
  int Nstrip = GammaRayTelDetector->GetNbOfTKRStrips();
  int Ntile = GammaRayTelDetector->GetNbOfTKRTiles();
  double sizexy = GammaRayTelDetector->GetTKRSizeXY();
  double sizez = GammaRayTelDetector->GetTKRSizeZ();
  int N = Nstrip*Ntile;      

  // Book1D histograms 
  //------------------

  // 1D histogram that store the energy deposition of the
  // particle in the last (number 0) TKR X-plane
  man->CreateH1("1","Edep in the last X plane (keV)", 100, 50, 200);

  // 1D histogram that store the hits distribution along the TKR X-planes
  man->CreateH1("2","Hits dist in TKR X planes",Nplane, 0, Nplane-1);

  // Book 2D histograms 
  //-------------------

  // 2D histogram that store the position (mm) of the hits (XZ projection)

  if (histo2DMode == "strip")
    {
      man->CreateH2("1","Tracker Hits XZ (strip,plane)", 
		    N, 0, N-1, 
		    2*Nplane, 0, Nplane-1); 
    }
  else
    {
      man->CreateH2("1","Tracker Hits XZ (x,z) in mm", 
		    int(sizexy/5), -sizexy/2, sizexy/2, 
		    int(sizez/5), -sizez/2, sizez/2);
    }  

  // 2D histogram that store the position (mm) of the hits (YZ projection)  
  if (histo2DMode == "strip")
    {
      man->CreateH2("2","Tracker Hits YZ (strip,plane)", 
		    N, 0, N-1, 
		    2*Nplane, 0, Nplane-1); 
    }
  else
    {
      man->CreateH2("2","Tracker Hits YZ (x,z) in mm", 
		    int(sizexy/5), -sizexy/2, sizexy/2, 
		    int(sizez/5), -sizez/2, sizez/2);
    }  
  
    
  // Book Ntuples (energy / plane/ x / y / z)
  //------------------------------------------  
  man->CreateNtuple("1","Track ntuple");
  man->CreateNtupleDColumn("energy");
  man->CreateNtupleDColumn("plane"); // can I use Int values? 
  man->CreateNtupleDColumn("x"); // can I use Int values?
  man->CreateNtupleDColumn("y"); // can I use Int values?
  man->CreateNtupleDColumn("z"); // can I use Int values?
  man->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   This member is called at the end of each run 
*/
void GammaRayTelAnalysis::EndOfRun() 
{
  //Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
}

/* This member is called at the end of every event */

void GammaRayTelAnalysis::EndOfEvent(G4int /* flag */ ) 
{;}








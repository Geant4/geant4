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
/// @file Analysis.cc
/// @brief Define histograms and ntuples

#include "G4AutoDelete.hh"
#include "G4SystemOfUnits.hh"
#include "Analysis.hh"

//Select format of output here
//Note: ntuple merging is supported only with Root format
#include "g4root.hh"
#include "G4RootAnalysisManager.hh"


G4ThreadLocal G4int Analysis::fincidentFlag = false;
G4ThreadLocal Analysis* the_analysis = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Analysis*
Analysis::GetAnalysis()
{
  if (!the_analysis)
    {
      the_analysis = new Analysis();
      G4AutoDelete::Register(the_analysis);
    }
  return the_analysis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Analysis::Analysis() :
  fUseNtuple(true),
  fMergeNtuple(true),
  fincident_x_hist(0), fincident_map(0), fdose_hist(0), fdose_map(0),
      fdose_prof(0), fdose_map_prof(0), fdose_map3d(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
Analysis::Book()
{
  G4cout << "Analysis::Book start, fUseNtuple: " << fUseNtuple << G4endl;

  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  mgr->SetVerboseLevel(1);
#ifdef G4MULTITHREADED
  // MT ntuple merging
  mgr->SetNtupleMerging(fMergeNtuple);
#endif

  fincident_x_hist = mgr->CreateH1("incident_x", "Incident X", 100, -5 * cm,
      5 * cm, "cm");
  fincident_map = mgr->CreateH2("incident_map", "Incident Map", 50, -5 * cm,
      5 * cm, 50, -5 * cm, 5 * cm, "cm", "cm");
  fdose_hist
      = mgr->CreateH1("dose", "Dose distribution", 500, 0, 50 * cm, "cm");
  fdose_map = mgr->CreateH2("dose_map", "Dose distribution", 500, 0, 50 * cm,
      200, -10 * cm, 10 * cm, "cm", "cm");
  fdose_map3d = mgr->CreateH3("dose_map_3d", "Dose distribution", 30, 0,
      50 * cm, 20, -10 * cm, 10 * cm, 20, -10 * cm, 10 * cm, "cm", "cm", "cm");
  fdose_prof = mgr->CreateP1("dose_prof", "Dose distribution", 300, 0, 30 * cm,
      0, 100 * MeV, "cm", "MeV");
  fdose_map_prof = mgr->CreateP2("dose_map_prof", "Dose distribution", 300, 0,
      30 * cm, 80, -4 * cm, 4 * cm, 0, 100 * MeV, "cm", "cm", "MeV");

  if ( fUseNtuple) {
    mgr->CreateNtuple("Dose", "Dose distribution");  // ntuple Id = 0
    mgr->CreateNtupleDColumn("pz");
    mgr->CreateNtupleDColumn("px");
    mgr->CreateNtupleDColumn("dose");
    mgr->FinishNtuple();
  }

  G4cout << "Analysis::Book finished " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Analysis::~Analysis()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
Analysis::OpenFile(const G4String& fname)
{
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  mgr->OpenFile(fname.c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
Analysis::Save()
{
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  mgr->Write();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
Analysis::Close(G4bool reset)
{
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  mgr->CloseFile(reset);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
Analysis::FillIncident(const G4ThreeVector& p)
{
  if (!fincidentFlag)
    {
      G4AnalysisManager* mgr = G4AnalysisManager::Instance();
      mgr->FillH2(fincident_map, p.x(), p.y());
      mgr->FillH1(fincident_x_hist, p.x());
      fincidentFlag = true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
Analysis::FillDose(const G4ThreeVector& p, G4double dedx)
{
  const G4double dxy = 10. * mm;
  if (std::abs(p.y()) < dxy)
    {
      G4AnalysisManager* mgr = G4AnalysisManager::Instance();
      const G4double Z0 = 25. * cm;

      mgr->FillH2(fdose_map, p.z() + Z0, p.x(), dedx / GeV);
      mgr->FillP2(fdose_map_prof, p.z() + Z0, p.x(), dedx);
      mgr->FillH3(fdose_map3d, p.z() + Z0, p.x(), p.y(), dedx / GeV);
      if (std::abs(p.x()) < dxy)
        {
          mgr->FillH1(fdose_hist, p.z() + Z0, dedx / GeV);
          mgr->FillP1(fdose_prof, p.z() + Z0, dedx);
        }

      if ( fUseNtuple ) { 
        // the same fill frequency as H2 "dose_map"
        mgr->FillNtupleDColumn(0, p.z() + Z0);
        mgr->FillNtupleDColumn(1, p.x());
        mgr->FillNtupleDColumn(2, dedx / GeV);
        mgr->AddNtupleRow();
      }
    }

}

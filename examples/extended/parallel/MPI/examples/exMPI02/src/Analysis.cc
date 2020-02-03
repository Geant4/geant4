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
/// @brief Define histograms

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "G4SystemOfUnits.hh"
#include "Analysis.hh"
#include "G4AutoLock.hh"

G4ThreadLocal G4int Analysis::fincidentFlag = false;
G4ThreadLocal Analysis* the_analysis = 0;

G4Mutex rootm = G4MUTEX_INITIALIZER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Analysis* Analysis::GetAnalysis()
{
   if ( ! the_analysis ) the_analysis = new Analysis;
   return the_analysis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Analysis::Analysis()
{
  G4AutoLock l(&rootm);
  // define histograms
  fincident_map = new TH2D("incident map", "Incident Distributuon",
                            50, -5., 5.,
                            50, -5., 5.);
  fincident_map-> GetXaxis()-> SetTitle("X (cm)");
  fincident_map-> GetYaxis()-> SetTitle("Y (cm)");
  fincident_map-> SetStats(0);

  fincident_x_hist = new TH1D("incident x", "Incident X", 100, -5., 5.);
  fincident_x_hist-> GetXaxis()-> SetTitle("X (cm)");
  fincident_x_hist-> SetFillColor(kRed);

  fdose_map = new TH2D("dose map", "Dose Distribution",
                       500, 0., 50.,
                       200, -10., 10.);
  fdose_map-> GetXaxis()-> SetTitle("Z (cm)");
  fdose_map-> GetYaxis()-> SetTitle("X (cm)");
  fdose_map-> SetStats(0);

  fdose_hist = new TH1D("dose", "Dose Distribution", 500, 0., 50.);
  fdose_hist-> GetXaxis()-> SetTitle("Z (cm)");
  fdose_hist-> GetYaxis()-> SetTitle("Dose (GeV)");
  fdose_hist-> SetFillColor(kBlue);
  fdose_hist-> SetStats(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Analysis::~Analysis()
{
  delete fincident_map;
  delete fincident_x_hist;
  delete fdose_map;
  delete fdose_hist;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Analysis::Update()
{
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Analysis::Clear()
{
  fincident_map-> Reset();
  fincident_x_hist-> Reset();
  fdose_map-> Reset();
  fdose_hist-> Reset();

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Analysis::Save(const G4String& fname)
{
  G4AutoLock l(&rootm);
  TFile* file = new TFile(fname.c_str(), "RECREATE", "Geant4 ROOT analysis");

  fincident_map-> Write();
  fincident_x_hist-> Write();
  fdose_map-> Write();
  fdose_hist-> Write();

  file-> Close();

  delete file;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Analysis::FillIncident(const G4ThreeVector& p)
{
  if ( ! fincidentFlag ) {
    fincident_map-> Fill(p.x()/cm, p.y()/cm);
    fincident_x_hist-> Fill(p.x()/cm);

    fincidentFlag = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Analysis::FillDose(const G4ThreeVector& p, G4double dedx)
{
  const G4double Z0 = 25.*cm;
  const G4double dxy = 10.*mm;

  if ( std::abs(p.y()) < dxy ) {
    fdose_map-> Fill((p.z()+Z0)/cm, p.x()/cm, dedx/GeV);

    if ( std::abs(p.x()) < dxy ) {
      fdose_hist-> Fill((p.z()+Z0)/cm, dedx/GeV);
    }
  }
}

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
/// @file Analysis.cc
/// @brief Define histograms

#include "Analysis.hh"
#include "G4SystemOfUnits.hh"

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

Analysis* Analysis::myanalysis = NULL;

// --------------------------------------------------------------------------
Analysis::Analysis()
{
  // ROOT style
  gROOT-> Reset();

  // define histograms
  incident_map = new TH2D("incident map", "Incident Distributuon",
                          50, -5., 5.,
                          50, -5., 5.);
  incident_map-> GetXaxis()-> SetTitle("X (cm)");
  incident_map-> GetYaxis()-> SetTitle("Y (cm)");
  incident_map-> SetStats(0);

  incident_x_hist = new TH1D("incident x", "Incident X", 100, -5., 5.);
  incident_x_hist-> GetXaxis()-> SetTitle("X (cm)");
  incident_x_hist-> SetFillColor(kRed);

  dose_map = new TH2D("dose map", "Dose Distribution", 
                      500, 0., 50.,
                      200, -10., 10.);
  dose_map-> GetXaxis()-> SetTitle("Z (cm)");
  dose_map-> GetYaxis()-> SetTitle("X (cm)");
  dose_map-> SetStats(0);

  dose_hist = new TH1D("dose", "Dose Distribution", 500, 0., 50.);
  dose_hist-> GetXaxis()-> SetTitle("Z (cm)");
  dose_hist-> GetYaxis()-> SetTitle("Dose (GeV)");
  dose_hist-> SetFillColor(kBlue);
  dose_hist-> SetStats(0);

}

// --------------------------------------------------------------------------
Analysis::~Analysis()
{
  delete incident_map;
  delete incident_x_hist;
  delete dose_map;
  delete dose_hist;

  myanalysis = NULL;
}

// --------------------------------------------------------------------------
Analysis* Analysis::GetAnalysis()
{
  if ( myanalysis == NULL ) {
    myanalysis = new Analysis();
  }

  return myanalysis;
}

// --------------------------------------------------------------------------
void Analysis::Update()
{
  return;
}

// --------------------------------------------------------------------------
void Analysis::Clear()
{
  incident_map-> Reset();
  incident_x_hist-> Reset();
  dose_map-> Reset();
  dose_hist-> Reset();

  return;
}

// --------------------------------------------------------------------------
void Analysis::Save(const G4String& fname)
{
  TFile* file = new TFile(fname.c_str(),
                          "RECREATE", "Geant4 ROOT analysis");
  incident_map-> Write();
  incident_x_hist-> Write();
  dose_map-> Write();
  dose_hist-> Write();

  file-> Close();

  delete file;

  return;
}

// --------------------------------------------------------------------------
void Analysis::FillIncident(const G4ThreeVector& p)
{
  if ( ! incidentFlag ) {
    incident_map-> Fill(p.x()/cm, p.y()/cm);
    incident_x_hist-> Fill(p.x()/cm);

    incidentFlag = true;
  }
}

// --------------------------------------------------------------------------
void Analysis::FillDose(const G4ThreeVector& p, G4double dedx)
{
  const G4double Z0 = 25.*cm;
  const G4double dxy = 10.*mm;

  if(std::abs(p.y()) < dxy ) {
    dose_map-> Fill((p.z()+Z0)/cm, p.x()/cm, dedx/GeV);
    
    if(std::abs(p.x()) < dxy) {
      dose_hist-> Fill((p.z()+Z0)/cm, dedx/GeV);
    }
  }
}


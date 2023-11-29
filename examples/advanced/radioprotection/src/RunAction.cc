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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it
// 

#include "RunAction.hh"
#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4Run.hh"

#ifdef ANALYSIS_USE
RunAction::RunAction(AnalysisManager* analysis, DetectorMessenger* detector)
{ 

  analysisMan = analysis;
  detectorMess = detector;

}
#else
RunAction::RunAction(DetectorMessenger* detector)
{
   detectorMess = detector;
}
#endif

RunAction::~RunAction()
{ }

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int run_number = aRun->GetRunID();
  G4cout << "### Run " << run_number << " start." << G4endl;

#ifdef ANALYSIS_USE
  G4bool additionalOutput;
  G4String detectorType = detectorMess -> GetDetectorType();
  if( detectorType == "DiamondTelescope" ) additionalOutput = true;
  else additionalOutput = false;
  
  // Create ROOT file, histograms and ntuple
  analysisMan -> book(additionalOutput);
#else
  G4cout << "Output to file disabled: no file will be generated. "
  << "To enable it, re-run cmake with the -DWITH_ANALYSIS_USE=ON flag." << G4endl;
#endif
  
  // check if there are pending changes in the geometry
  if( detectorMess -> AreTherePendingChanges() )
  {
    G4cout << "WARNING: pending changes to the geometry found during BeginOfRunAction. " 
    << "These will be ignored during the current run. "
    << "Use /geometrySetup/applyChanges if you want to apply them in your next run." << G4endl;
  }
  
  // print the current geometry setup
  G4cout << "Simulating the " << detectorMess->GetDetectorType() << " detector "
  << "with dimensions: " << detectorMess->GetDetectorSizeWidth()/um << " um width, "
  << detectorMess->GetDetectorSizeThickness()/um << " um thickness." << G4endl;
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "Number of events = " << aRun->GetNumberOfEvent() << G4endl;

#ifdef ANALYSIS_USE
// Close the output ROOT file with the results
   analysisMan -> finish(); 
#endif
}


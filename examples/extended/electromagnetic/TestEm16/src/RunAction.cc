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
// $Id: RunAction.cc,v 1.6 2006/06/29 16:48:02 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"

#ifdef G4ANALYSIS_USE
 #include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
:af(0), tree(0)
{
  for (G4int j=0; j<n_histos; j++) histo[j] = 0;
  n_gam_sync = 0;
  e_gam_sync = 0;
  e_gam_sync2 = 0;
  e_gam_sync_max =0;
  lam_gam_sync = 0;

#ifdef G4ANALYSIS_USE
 // Creating the analysis factory
 af = AIDA_createAnalysisFactory();

 if (af) {
   // Creating the tree factory
   AIDA::ITreeFactory* tf = af->createTreeFactory();

   // Creating a tree mapped to an hbook file.
   G4bool readOnly  = false;
   G4bool createNew = true;
   G4String options = "--noErrors uncompress";
   tree = tf->create("testem16.hbook","hbook",readOnly,createNew,options);
   //tree = tf->create("testem16.root", "root",readOnly,createNew,options);
   //tree = tf->create("testem16.xml" ,"xml"  ,readOnly,createNew,options);
   delete tf;

   if (tree) {
    // Creating a histogram factory
    AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);

    // Creating histograms
    const G4double Ecr=66.5025;
    histo[0] = hf->createHistogram1D("1","SynRad Energy in keV",100, 0 ,5.*Ecr);
    histo[1] = hf->createHistogram1D("2","SynRad Power  in keV",100, 0 ,5.*Ecr);
    histo[2] = hf->createHistogram1D("3","Path Length in m",100, 0, 1.6);

     delete hf;
     G4cout << "\n----> Histogram tree is opened" << G4endl;
   }
 }
  G4cout << "G4ANALYSIS_USE was     set, there will be    AIDA histos" << '\n';
#else
  G4cout << "G4ANALYSIS_USE was not set, there will be no AIDA histos" << '\n';
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{

#ifdef G4ANALYSIS_USE
  bool debug=true;
  bool Commit_Ok=tree->commit();       // Writing the histograms to the file
  if (Commit_Ok)
    { if(debug) G4cout << "tree->commit() ok. Writing of histogram file done" 
                       << '\n';
    }
  else if(!Commit_Ok)
    { G4cout << "tree->commit() not successful, no histogram file written" 
                      << '\n';
    }
  tree->close();        // and closing the tree (and the file)

  delete tree;
  delete af;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  if(n_gam_sync>0)
  {
    G4double Emean = e_gam_sync/n_gam_sync;
    G4double E_rms = std::sqrt(e_gam_sync2/n_gam_sync - Emean*Emean);
    G4cout
    << "Summary for synchrotron radiation :" << '\n' << std::setprecision(4)
    << "  Number of photons = " << n_gam_sync << '\n'
    << "  Emean             = " << Emean/keV << " +/- "
    << E_rms/(keV * std::sqrt((G4double) n_gam_sync)) << " keV" << '\n'
    << "  E_rms             = " << G4BestUnit(E_rms,"Energy") << '\n'
    << "  Energy Max / Mean = " << e_gam_sync_max / Emean << '\n'
    << "  MeanFreePath      = " << G4BestUnit(lam_gam_sync/n_gam_sync,"Length")
    << G4endl;
  }
  
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

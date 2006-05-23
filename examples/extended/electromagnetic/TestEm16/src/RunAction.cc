//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: RunAction.cc,v 1.3 2006-05-23 16:13:01 hbu Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
   //tree = tf->create("testem16.xml" ,"xml"  ,readOnly,createNew,options); // standard AIDA
   delete tf;

   if (tree) {
     // Creating a histogram factory
     AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);

     // Creating histograms
     static const G4double Ecr=66.5025;
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
  if(Commit_Ok)
  { if(debug) G4cout << "tree->commit() ok. Writing of histogram file done" << '\n';
  }
  else if(!Commit_Ok)
  { G4cout << "tree->commit() not successful, no histogram file written" << '\n';
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
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();

  if(n_gam_sync>0)
  {
    G4double Emean = e_gam_sync/n_gam_sync;
    G4double E_rms = sqrt(e_gam_sync2/n_gam_sync - Emean*Emean);
    G4cout
    << "Summary for synchrotron radiation :" << '\n' << std::setprecision(4)
    << "  Number of photons = " << n_gam_sync << '\n'
    << "  Emean             = " << Emean/keV << " +/- "
    << E_rms/(keV * sqrt((G4double) n_gam_sync)) << " keV" << '\n'
    << "  E_rms             = " << G4BestUnit(E_rms,"Energy") << '\n'
    << "  Energy Max / Mean = " << e_gam_sync_max / Emean << '\n'
    << "  MeanFreePath      = " << G4BestUnit(lam_gam_sync/n_gam_sync,"Length") << '\n';
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

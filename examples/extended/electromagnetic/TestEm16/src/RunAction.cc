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
// $Id: RunAction.cc,v 1.8 2007-01-18 10:24:35 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(HistoManager* histo)
:histoManager(histo)
{
  n_gam_sync = 0;
  e_gam_sync = 0;
  e_gam_sync2 = 0;
  e_gam_sync_max =0;
  lam_gam_sync = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();

  //book histograms      
  histoManager->book();
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

  //save histograms      
  histoManager->save();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

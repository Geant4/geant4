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
/// \file hadronic/Hadr00/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.hh,v 1.1 2008-07-07 16:37:26 vnivanch Exp $
//
// -------------------------------------------------------------
//  
//    GEANT4 class file
//    RunAction
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "HistoManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(const DetectorConstruction* det) 
  : fDetector(det)
{
  fHisto = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fHisto;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4bool show = true;
  //#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread() == true) { show = false; }
  //#endif
  //  G4cout << "###### show= " << show << G4endl;

  if(show) {
    G4int id = aRun->GetRunID();
    G4cout << "### Run " << id << " start" << G4endl;
    fHisto->BeginOfRun();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  G4bool show = true;
  //#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread() == true) { show = false; }
  //#endif

  if(show) {
    G4cout << "RunAction: End of run actions are started" << G4endl;
    if(fDetector) {
      fHisto->SetTargetMaterial(fDetector->GetTargetMaterial());
    }
    fHisto->EndOfRun();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

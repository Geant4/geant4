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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 78723 2014-01-20 10:32:17Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* event, DetectorConstruction* detector)
:G4UserSteppingAction(),
<<<<<<< HEAD
 fEventAction(event)
{ }
=======
 fEventAction(event),
 fDetectorConstruction(detector)
{}
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
 G4double edep = aStep->GetTotalEnergyDeposit();
 if (edep <= 0.) return;
 
 //total energy deposit in cytoplasm or nucleus
 //
 if (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()
     ==fDetectorConstruction->GetCytoLogicalVolume()) fEventAction->AddCytoEdep(edep);
       
 if (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()
     ==fDetectorConstruction->GetNuclLogicalVolume()) fEventAction->AddNuclEdep(edep);     

 //G4cout << edep << G4endl;

 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
     
 //step size of primary particle or charged secondaries
 //
 G4double steplen = aStep->GetStepLength();
 const G4Track* track = aStep->GetTrack();
 if      (track->GetTrackID() == 1) analysisManager->FillH1(8, steplen);
 else if (track->GetDefinition()->GetPDGCharge() != 0.)
                                    analysisManager->FillH1(9, steplen); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction(),
 fEdep1(0.), fEdep2(0.), fWeight1(0.), fWeight2(0.),
 fTime0(-1*s)
{ } 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep1 = fEdep2 = fWeight1 = fWeight2 = 0.;
  fTime0 = -1*s;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdep(G4int iVol, G4double edep,
                                      G4double time, G4double weight)
{
  // initialize t0
  if (fTime0 < 0.) fTime0 = time;
  
  // out of time window ?
  const G4double TimeWindow (1*microsecond);
  if (std::fabs(time - fTime0) > TimeWindow) return;
  
  if (iVol == 1) { fEdep1 += edep; fWeight1 += edep*weight;}
  if (iVol == 2) { fEdep2 += edep; fWeight2 += edep*weight;}  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
 G4double Etot = fEdep1 + fEdep2;
 G4double Wtot = (fWeight1 + fWeight2)/Etot;
 
 // pulse height in target
 //
 if (fEdep1 > 0.) {
   fWeight1 /= fEdep1;
   analysisManager->FillH1(0, fEdep1, fWeight1);   
 }
 
 // pulse height in detector
 //   
 if (fEdep2 > 0.) {
   fWeight2 /= fEdep2;
   analysisManager->FillH1(1, fEdep2, fWeight2);
 }
   
 // total
 //
 analysisManager->FillH1(2, Etot, Wtot);
 
 // threshold in target and detector        
 const G4double Threshold1(10*keV), Threshold2(10*keV);
  
 //coincidence, anti-coincidences 
 //  
 G4bool coincidence       = ((fEdep1 >= Threshold1) && (fEdep2 >= Threshold2));
 G4bool anti_coincidence1 = ((fEdep1 >= Threshold1) && (fEdep2 <  Threshold2));
 G4bool anti_coincidence2 = ((fEdep1 <  Threshold1) && (fEdep2 >= Threshold2)); 

 if (coincidence)       analysisManager->FillH1(3, fEdep2, fWeight2);
 if (anti_coincidence1) analysisManager->FillH1(4, fEdep1, fWeight1);
 if (anti_coincidence2) analysisManager->FillH1(5, fEdep2, fWeight2); 

 // pass energies to Run
 //  
 Run* run = static_cast<Run*>(
            G4RunManager::GetRunManager()->GetNonConstCurrentRun());
             
 run->AddEdep (fEdep1, fEdep2);             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



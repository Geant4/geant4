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
/// \file electromagnetic/TestEm5/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* DET, EventAction* EA)
:G4UserTrackingAction(),fDetector(DET), fEventAction(EA)
{ 
  fXstartAbs = fXendAbs = fPrimaryCharge = fDirX = 0.0;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack )
{
  // few initialisations
  //
  if (aTrack->GetTrackID() == 1) {
    fXstartAbs = fDetector->GetxstartAbs();
    fXendAbs   = fDetector->GetxendAbs();
    fPrimaryCharge = aTrack->GetDefinition()->GetPDGCharge();
    fDirX = aTrack->GetMomentumDirection().x();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 

  G4ThreeVector position = aTrack->GetPosition();
  G4ThreeVector vertex   = aTrack->GetVertexPosition();  
  G4double charge        = aTrack->GetDefinition()->GetPDGCharge();
  G4double energy        = aTrack->GetKineticEnergy();

  G4bool transmit = (((fDirX >= 0.0 && position.x() >= fXendAbs) ||
                      (fDirX < 0.0 && position.x() <= fXstartAbs))
                      && energy > 0.0);
  G4bool reflect  = (((fDirX >= 0.0 && position.x() <= fXstartAbs) ||
                      (fDirX < 0.0 && position.x() <= fXendAbs))
                      && energy > 0.0);
  G4bool notabsor = (transmit || reflect);

  //transmitted + reflected particles counter
  //
  G4int flag = 0;
  if (charge == fPrimaryCharge)  flag = 1;
  if (aTrack->GetTrackID() == 1) flag = 2;
  if (transmit) fEventAction->SetTransmitFlag(flag);
  if (reflect)  fEventAction->SetReflectFlag(flag);
      
  //
  //histograms
  //
  G4bool charged = (charge != 0.);
  G4bool neutral = !charged;

  //energy spectrum at exit
  //zero energy charged particle are not taken into account
  // 
  G4int id = 0;  
       if (transmit && charged) id = 10;
  else if (transmit && neutral) id = 20;
  else if (reflect  && charged) id = 30;
  else if (reflect  && neutral) id = 40;

  if (id>0) { analysisManager->FillH1(id, energy); }

  //energy leakage
  //
  if (notabsor) {
    G4int trackID = aTrack->GetTrackID();
    G4int index = 0; if (trackID > 1) index = 1;    //primary=0, secondaries=1
    G4double eleak = aTrack->GetKineticEnergy();
    if ((aTrack->GetDefinition() == G4Positron::Positron()) && (trackID > 1))
      eleak += 2*electron_mass_c2;
    run->AddEnergyLeak(eleak,index);  
  }

  //space angle distribution at exit : dN/dOmega
  //
  G4ThreeVector direction = aTrack->GetMomentumDirection();    
  id = 0;   
       if (transmit && charged) id = 12;
  else if (transmit && neutral) id = 22;
  else if (reflect  && charged) id = 32;
  else if (reflect  && neutral) id = 42;

  if (id>0) {
    G4double theta  = std::acos(direction.x());
    if (theta > 0.0) {
      G4double dteta  = analysisManager->GetH1Width(id);
      G4double unit   = analysisManager->GetH1Unit(id);    
      G4double weight = (unit*unit)/(twopi*std::sin(theta)*dteta);
      /*
      G4cout << "theta, dteta, unit, weight: " 
             << theta << "  "   
             << dteta << "  "   
             << unit << "  "   
             << weight << G4endl;   
      */
      analysisManager->FillH1(id,theta,weight);
    }
  }
  
  //energy fluence at exit : dE(MeV)/dOmega
  //
  id = 0;  
       if (transmit && charged) id = 11;
  else if (reflect  && charged) id = 31;
  else if (transmit && neutral) id = 21;
  else if (reflect  && neutral) id = 41;

  if (id>0) {
    G4double theta  = std::acos(direction.x());
    if (theta > 0.0) {
      G4double dteta  = analysisManager->GetH1Width(id);
      G4double unit   = analysisManager->GetH1Unit(id);    
      G4double weight = (unit*unit)/(twopi*std::sin(theta)*dteta);
      weight *= (aTrack->GetKineticEnergy()/MeV); 
      analysisManager->FillH1(id,theta,weight);
    } 
  }
  
  //projected angles distribution at exit
  //
  id = 0;   
       if (transmit && charged) id = 13;
  else if (transmit && neutral) id = 23;
  else if (reflect  && charged) id = 33;
  else if (reflect  && neutral) id = 43;

  if (id>0) {
    if (direction.x() != 0.0) {
      G4double tet = std::atan(direction.y()/std::fabs(direction.x()));
      analysisManager->FillH1(id,tet);
      if (transmit && (flag == 2)) run->AddMscProjTheta(tet);

      tet = std::atan(direction.z()/std::fabs(direction.x()));
      analysisManager->FillH1(id,tet);
      if (transmit && (flag == 2)) run->AddMscProjTheta(tet);
    }
  }

  //projected position and radius at exit
  //
  if (transmit && energy > 0.0) {
    G4double y = position.y(), z = position.z();
    G4double r = std::sqrt(y*y + z*z);
    analysisManager->FillH1(14,  y);
    analysisManager->FillH1(14,  z);
    analysisManager->FillH1(15,  r);
  }
  
  //x-vertex of charged secondaries
  //
  if ((aTrack->GetParentID() == 1) && charged && energy > 0.0) {
    G4double xVertex = (aTrack->GetVertexPosition()).x();
    analysisManager->FillH1(6, xVertex);
    if (notabsor) analysisManager->FillH1(7, xVertex); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


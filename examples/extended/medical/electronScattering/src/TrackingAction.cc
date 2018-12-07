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
/// \file medical/electronScattering/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4Track.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* DET)
:fDetector(DET)
{
 fZend = 0.5*(fDetector->GetThicknessWorld()); 
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track*)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  G4double charge = track->GetDefinition()->GetPDGCharge();
  G4ThreeVector position = track->GetPosition();
  G4ThreeVector direction = track->GetMomentumDirection();    

  if (charge == 0.0)       return;
  if (position.z() < fZend) return;
  if (direction.z() <= 0.) return;
    
  G4double rmin, dr, ds;
  G4int ih = 1;
    
  //projected angle at exit
  
  G4double ux = direction.x(), uy = direction.y(), uz = direction.z();
  G4double thetax = std::atan(ux/uz);
  G4double thetay = std::atan(uy/uz);
  analysisManager->FillH1(ih, thetax);
  analysisManager->FillH1(ih, thetay);
      
  //dN/dS at exit
  
  G4double x = position.x(), y = position.y();
  G4double r = std::sqrt(x*x + y*y);
  ih = 2;
  dr = analysisManager->GetH1Width(ih);
  rmin = ((int)(r/dr))*dr;
  ds = twopi*(rmin + 0.5*dr)*dr;        
  analysisManager->FillH1(ih, r, 1/ds);  
      
  //d(N/cost)/dS at exit
  
  ih = 3;
  dr = analysisManager->GetH1Width(ih);
  rmin = ((int)(r/dr))*dr;
  ds = twopi*(rmin + 0.5*dr)*dr;        
  analysisManager->FillH1(ih, r, 1/(uz*ds));
  
  //vector of d(N/cost)/dS at exit
  
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  run->SumFluence(r, 1/uz);
  
  //space angle at exit : dN/dOmega
  
  ih = 5;
  G4double theta = std::acos(uz);
  if (theta > 0.) {
     G4double dtheta = analysisManager->GetH1Width(ih);
     G4double unit   = analysisManager->GetH1Unit(ih);
     if (dtheta > 0.) {
        G4double weight = unit*unit/(twopi*std::sin(theta)*dtheta);
        analysisManager->FillH1(ih, theta, weight);
     }
  }
  
  //measured space angle at exit : dN/dOmega_meas
  //
  ih = 6;  
  const G4double dist = fDetector->GetZdist_foil_detector();
  G4double thetam = std::atan(r/dist);
  if (thetam > 0.) {
     G4double dtheta = analysisManager->GetH1Width(ih);
     G4double unit   = analysisManager->GetH1Unit(ih);
     if (dtheta > 0.) {
        G4double weight = unit*unit/(twopi*std::sin(thetam)*dtheta);
        analysisManager->FillH1(ih, thetam, weight);
     }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


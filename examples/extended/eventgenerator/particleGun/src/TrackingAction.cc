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
/// \file eventgenerator/particleGun/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction2.hh"
#include "PrimaryGeneratorAction3.hh"
#include "PrimaryGeneratorAction4.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(PrimaryGeneratorAction* prim)
:G4UserTrackingAction(), fPrimary(prim)
{ 
 // parameters for generator action #3
  fNewUz = fPrimary->GetAction3()->GetNewUz();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{ 
 G4int selectedGeneratorAction = fPrimary->GetSelectedAction();
 G4AnalysisManager* analysis = G4AnalysisManager::Instance(); 
 G4int id = 0;

 if(selectedGeneratorAction==2)
 {
  //energy spectrum
  //
  id = 1;   
  analysis->FillH1(id, track->GetKineticEnergy());
 }
 
 else if(selectedGeneratorAction==0)
 {
  //particle direction: cos(alpha)
  //
  id = 5;
  G4ThreeVector um = track->GetMomentumDirection();
  G4double cosalpha  = um.z();
  analysis->FillH1(id, cosalpha);
      
  //particle direction: psi
  //
  id = 6;
  G4double psi = std::atan2(um.y(), um.x());
  if (psi < 0.) psi += twopi;    
  analysis->FillH1(id, psi);  
 }
 else if(selectedGeneratorAction==3)
 {
  //particle direction in local frame: cos(alpha)
  //
  id = 5;
  G4ThreeVector um = track->GetMomentumDirection();
  G4double cosalpha  = fNewUz*um;
  analysis->FillH1(id, cosalpha);
      
  //particle direction in local frame: psi
  //
  id = 6;
  // complete local frame
  G4ThreeVector u1(1.,0.,0.);  u1.rotateUz(fNewUz);
  G4ThreeVector u2(0.,1.,0.);  u2.rotateUz(fNewUz);  
  //  
  G4double psi = std::atan2(u2*um, u1*um);
  if (psi < 0.) psi += twopi;    
  analysis->FillH1(id, psi);  
 }

 else if(selectedGeneratorAction==4)
 {
  G4ThreeVector vertex   = track->GetVertexPosition();
  G4double r = vertex.mag();
  if (r <= 0.0) return;
  //local frame : new uz = ur
  G4ThreeVector ur   = vertex/r;

  //vertex position. radial distribution dN/dv = f(r)
  //
  id = 2;
  G4double dr = analysis->GetH1Width(id);
  G4double dv = 2*twopi*r*r*dr;
  if (dv > 0.) analysis->FillH1(id, r, 1./dv);

  //vertex position: cos(theta)
  //
  id = 3;
  G4double costheta  = ur.z();
  analysis->FillH1(id, costheta);

  //vertex position: phi
  //
  id = 4;
  G4double phi  = std::atan2(ur.y(), ur.x());
  if (phi < 0.) phi += twopi;
  analysis->FillH1(id, phi);

  //particle direction in local frame: cos(alpha)
  //
  id = 5;
  G4ThreeVector um = track->GetMomentumDirection();
  G4double cosalpha  = ur*um;
  analysis->FillH1(id, cosalpha);

  //particle direction in local frame: psi
  //
  id = 6;
  // complete local frame
  G4ThreeVector u1(1.,0.,0.);  u1.rotateUz(ur);
  G4ThreeVector u2(0.,1.,0.);  u2.rotateUz(ur);
  //
  G4double psi = std::atan2(u2*um, u1*um);
  if (psi < 0.) psi += twopi;
  analysis->FillH1(id, psi);
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track*)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


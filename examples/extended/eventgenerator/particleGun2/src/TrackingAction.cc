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
//
// $Id: TrackingAction.cc,v 1.1 2010-05-16 21:32:31 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "TrackingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(PrimaryGeneratorAction* prim,HistoManager* histo)
:primary(prim),histoManager(histo)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  //energy spectrum
  //
  G4int id = 1;   
  histoManager->FillHisto(id, track->GetKineticEnergy());

  //local frame : new uz = ur
  G4ThreeVector newUz = primary->GetNewUz();
    
  //momentum direction : angular distr dN/dOmega = f(alpha)
  //
  id = 2;
  G4ThreeVector um = track->GetMomentumDirection();
  G4double alpha  = std::acos(newUz*um);
  G4double dalfa  = histoManager->GetBinWidth(id);
  G4double dOmega = twopi*std::sin(alpha)*dalfa;     
  if (dOmega > 0.) histoManager->FillHisto(id, alpha, 1./dOmega);
    
  //momentum direction : angular distr dN/dOmega = f(psi)
  //
  id = 3;
  // complete local frame
  G4ThreeVector u1(1.,0.,0.);  u1.rotateUz(newUz);
  G4ThreeVector u2(0.,1.,0.);  u2.rotateUz(newUz);  
  //  
  G4double psi = std::atan2(u2*um, u1*um);
  if (psi < 0.) psi += twopi;    
  G4double dpsi  = histoManager->GetBinWidth(id);
  G4double alphaMax = primary->GetAlphaMax();    
  dOmega = (1. - std::cos(alphaMax))*dpsi;     
  if (dOmega > 0.) histoManager->FillHisto(id, psi, 1./dOmega);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track*)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


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
// $Id: TrackingAction.cc,v 1.11 2006/06/29 16:56:29 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* DET,RunAction* RA,
                               EventAction* EA, HistoManager* HM)
:detector(DET), runaction(RA), eventaction(EA), histoManager(HM)
{ }
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack )
{
  // few initialisations
  //
  if (aTrack->GetTrackID() == 1) {
    xstartAbs = detector->GetxstartAbs();
    xendAbs   = detector->GetxendAbs();
    primaryCharge = aTrack->GetDefinition()->GetPDGCharge();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4ThreeVector position = aTrack->GetPosition();
  G4double charge    = aTrack->GetDefinition()->GetPDGCharge();

  G4bool transmit = (position.x() >= xendAbs);
  G4bool reflect  = (position.x() <= xstartAbs);

  //transmitted + reflected particles counter
  //
  G4int flag = 0;
  if (charge == primaryCharge)   flag = 1;
  if (aTrack->GetTrackID() == 1) flag = 2;
  if (transmit) eventaction->SetTransmitFlag(flag);
  if (reflect)  eventaction->SetReflectFlag(flag);

  //
  //histograms
  //
  G4int id = 0;
  G4bool charged  = (charge != 0.);
  G4bool neutral = !charged;

  //energy spectrum at exit
  //
       if (transmit && charged) id =  4;
  else if (transmit && neutral) id =  8;
  else if (reflect  && charged) id = 11;
  else if (reflect  && neutral) id = 14;

  histoManager->FillHisto(id, aTrack->GetKineticEnergy());

  //space angle distribution at exit
  //
       if (transmit && charged) id =  5;
  else if (transmit && neutral) id =  9;
  else if (reflect  && charged) id = 12;
  else if (reflect  && neutral) id = 15;

  G4ThreeVector direction = aTrack->GetMomentumDirection();
  if (histoManager->HistoExist(id)) {
    G4double theta  = std::acos(direction.x());
    G4double dteta  = histoManager->GetBinWidth(id);
    G4double weight = 1./(twopi*std::sin(theta)*dteta);
    G4double unit   = histoManager->GetHistoUnit(id);
    histoManager->FillHisto(id,theta,weight*unit*unit);
  }

  //projected angles distribution at exit
  //
       if (transmit && charged) id =  6;
  else if (transmit && neutral) id = 10;
  else if (reflect  && charged) id = 13;
  else if (reflect  && neutral) id = 16;

  if(id>0) {
    G4double tet = std::atan(direction.y()/std::fabs(direction.x()));
    histoManager->FillHisto(id,tet);
    if (transmit && (flag == 2)) runaction->AddMscProjTheta(tet);

    tet = std::atan(direction.z()/std::fabs(direction.x()));
    histoManager->FillHisto(id,tet);
    if (transmit && (flag == 2)) runaction->AddMscProjTheta(tet);
  }

  //projected position at exit
  //
  if (transmit && charged) {
    histoManager->FillHisto(7, position.y());
    histoManager->FillHisto(7, position.z());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


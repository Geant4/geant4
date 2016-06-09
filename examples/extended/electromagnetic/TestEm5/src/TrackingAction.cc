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

//
// $Id: TrackingAction.cc,v 1.2 2003/08/27 17:18:16 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"
 
#ifdef G4ANALYSIS_USE
 #include "AIDA/IHistogram1D.h"
#endif

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
    worldLimit = 0.5*(detector->GetWorldSizeX());
    primaryCharge = aTrack->GetDefinition()->GetPDGCharge();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4ThreeVector position = aTrack->GetPosition();
  G4double charge    = aTrack->GetDefinition()->GetPDGCharge();

  G4bool transmit = (position.x() >=  worldLimit);
  G4bool reflect  = (position.x() <= -worldLimit);

  //transmitted + reflected particles counter
  //
  G4int flag = 0;
  if (charge == primaryCharge)   flag = 1;
  if (aTrack->GetTrackID() == 1) flag = 2;
  if (transmit) eventaction->SetTransmitFlag(flag);
  if (reflect)  eventaction->SetReflectFlag(flag);
  
#ifdef G4ANALYSIS_USE
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
  if (histoManager->GetHisto(id)) {
    G4double energy = aTrack->GetKineticEnergy();
    G4double unit   = histoManager->GetHistoUnit(id);
    histoManager->GetHisto(id)->fill(energy/unit);
  }
  
  //space angle distribution at exit
  //
       if (transmit && charged) id =  5;
  else if (transmit && neutral) id =  9;
  else if (reflect  && charged) id = 12;
  else if (reflect  && neutral) id = 15;
  if (histoManager->GetHisto(id)) {
    G4ThreeVector direction = aTrack->GetMomentumDirection();
    G4double theta  = acos(abs(direction.x()));
    G4double dteta  = histoManager->GetBinWidth(id);
    G4double weight = 1./(2*pi*sin(theta)*dteta);  
    G4double unit   = histoManager->GetHistoUnit(id);
    histoManager->GetHisto(id)->fill(theta/unit,weight*unit*unit); 
  }
    
  //projected angle distribution at exit
  //
       if (transmit && charged) id =  6;
  else if (transmit && neutral) id = 10;
  else if (reflect  && charged) id = 13;
  else if (reflect  && neutral) id = 16;
  if (histoManager->GetHisto(id)) {
    G4ThreeVector momentum = aTrack->GetMomentum();  
    G4double unit = histoManager->GetHistoUnit(id);
    histoManager->GetHisto(id)->fill(atan(momentum.y()/abs(momentum.x()))/unit);
    histoManager->GetHisto(id)->fill(atan(momentum.z()/abs(momentum.x()))/unit);
  }
            
  //radial dispersion at exit
  //
  if (transmit && charged) id =  7;
  if (histoManager->GetHisto(id)) {
    G4double radius = sqrt(position.y()*position.y()+position.z()*position.z());
    G4double unit   = histoManager->GetHistoUnit(id);
    histoManager->GetHisto(id)->fill(radius/unit);
  }  
#endif    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


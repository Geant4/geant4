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
// $Id: TrackingAction.cc,v 1.2 2010-05-22 21:21:52 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* DET,RunAction* RA,
                               HistoManager* HM)
:detector(DET), runaction(RA), histoManager(HM)
{
 Zend = 0.5*(detector->GetThicknessWorld()); 
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track*)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  G4double charge = track->GetDefinition()->GetPDGCharge();
  G4ThreeVector position = track->GetPosition();
  G4ThreeVector direction = track->GetMomentumDirection();    

  if (charge == 0.0)       return;
  if (position.z() < Zend) return;
  if (direction.z() <= 0.) return;
    
  G4double rmin, dr, ds;
  G4int ih = 1;
    
  //projected angle at exit
  //
  G4double ux = direction.x(), uy = direction.y(), uz = direction.z();
  if (histoManager->HistoExist(ih)) {
    G4double thetax = std::atan(ux/uz);
    G4double thetay = std::atan(uy/uz);
    histoManager->FillHisto(ih, thetax);
    histoManager->FillHisto(ih, thetay);
  }
      
  //dN/dS at exit
  //
  G4double x = position.x(), y = position.y();
  G4double r = std::sqrt(x*x + y*y);
  ih = 2;
  if (histoManager->HistoExist(ih)) {
    dr = histoManager->GetBinWidth(ih);
    rmin = ((int)(r/dr))*dr;
    ds = twopi*(rmin + 0.5*dr)*dr;        
    histoManager->FillHisto(ih, r, 1/ds);  
  }
      
  //d(N/cost)/dS at exit
  //
  ih = 3;
  if (histoManager->HistoExist(ih)) {
    dr = histoManager->GetBinWidth(ih);
    rmin = ((int)(r/dr))*dr;
    ds = twopi*(rmin + 0.5*dr)*dr;        
    histoManager->FillHisto(ih, r, 1/(uz*ds));
  }
  
  //vector of d(N/cost)/dS at exit
  //
  runaction->SumFluence(r, 1/uz);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


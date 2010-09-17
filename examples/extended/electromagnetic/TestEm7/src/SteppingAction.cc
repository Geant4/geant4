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
// $Id: SteppingAction.cc,v 1.15 2010-09-17 18:45:43 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, HistoManager* histo,
                                RunAction* RuAct)
:detector(det), histoManager(histo), runAction(RuAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.) return;
  
  G4double niel = step->GetNonIonizingEnergyDeposit();

  runAction->FillEdep(edep, niel);

  if (step->GetTrack()->GetTrackID() == 1) runAction->AddPrimaryStep();  
 
  //Bragg curve
  //	
  G4StepPoint* prePoint  = step->GetPreStepPoint();
  G4StepPoint* postPoint = step->GetPostStepPoint();
   
  G4double x1 = prePoint->GetPosition().x();
  G4double x2 = postPoint->GetPosition().x();  
  G4double x  = x1 + G4UniformRand()*(x2-x1) + 0.5*(detector->GetAbsorSizeX());
  histoManager->FillHisto(1, x, edep);
  histoManager->FillHisto(2, x, edep);

  //fill tallies
  //
  G4TouchableHandle touchable = prePoint->GetTouchableHandle();
  G4LogicalVolume* lVolume = touchable->GetVolume()->GetLogicalVolume();
  if (lVolume == detector->GetLogicalTally())
    runAction->FillTallyEdep(touchable->GetCopyNumber(), edep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



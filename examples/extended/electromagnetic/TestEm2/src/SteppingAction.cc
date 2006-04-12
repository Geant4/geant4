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
// $Id: SteppingAction.cc,v 1.5 2006-04-12 15:25:38 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"

#include "G4SteppingManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* run)
:detector(det),Run(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
 // energy deposit
 //
 G4double dEStep = step->GetTotalEnergyDeposit();
 if (dEStep > 0.) {
   G4ThreeVector prePoint  = step->GetPreStepPoint()->GetPosition();
   G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
   G4ThreeVector position  = prePoint + G4UniformRand()*(postPoint - prePoint);
   G4double x = position.x(), y = position.y(), z = position.z();
   G4double radius = std::sqrt(x*x + y*y);
   G4double offset = 0.5*detector->GetfullLength();
   G4int SlideNb = int((z + offset)/detector->GetdLlength());
   G4int RingNb  = int(radius/detector->GetdRlength());        
   Run->fillPerStep(dEStep,SlideNb,RingNb);
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
